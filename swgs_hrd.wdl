version 1.0

import "BamMetrics/bammetrics.wdl" as bammetrics
import "structs.wdl" as structs
import "QC/QC.wdl" as qc
import "tasks/biowdl.wdl" as biowdl
import "tasks/bwa.wdl" as bwa
import "tasks/freec.wdl" as freec
import "tasks/gem2.wdl" as gem2
import "tasks/multiqc.wdl" as multiqc
import "tasks/picard.wdl" as picard
import "tasks/samtools.wdl" as samtools
import "tasks/shallowHRD.wdl" as shallowHRD


workflow swgs_hrd {
    input {
        File sampleConfigFile
        BwaIndex bwaIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict

        String adapterForward = "AGATCGGAAGAG"
        String adapterReverse = "AGATCGGAAGAG"
        String platform = "illumina"
        Int readLength = 150

        File? mappability
        Array[File]+? chrFiles
    }

    meta {
        allowNestedInputs: true
    }

    # mappbility
    if (! defined(mappability)) {
        call gem2.Index as gemIndex {
            input:
                fasta = referenceFasta
        }

        call gem2.Mappability as gemMappability {
            input:
                gemIndex = gemIndex.gemIndex,
                readLength = readLength
        }
    }

    if (! defined(chrFiles)) {
        call SplitFasta as splitFasta {
            input:
                fasta = referenceFasta,
                outDir = "./chrFiles"
        }
    }

    call biowdl.InputConverter as convertSampleConfig {
        input:
            samplesheet = sampleConfigFile
    }

    SampleConfig sampleConfig = read_json(convertSampleConfig.json)

    scatter (sample in sampleConfig.samples) {
        String sampleDir = "./~{sample.id}"
        scatter (readgroup in sample.readgroups) {
            String readgroupDir = "~{sampleDir}/~{readgroup.id}"

            call qc.QC as qualityControl {
                input:
                    outputDir = readgroupDir,
                    read1 = readgroup.R1,
                    read2 = readgroup.R2,
                    adapterForward = adapterForward,
                    adapterReverse = adapterReverse
            }

            call bwa.Mem as bwaMem {
                input:
                    read1 = qualityControl.qcRead1,
                    read2 = qualityControl.qcRead2,
                    outputPrefix = readgroupDir + "/" + sample.id + "-" + readgroup.lib_id + "-" + readgroup.id,
                    readgroup = "@RG\\tID:~{sample.id}-~{readgroup.lib_id}-~{readgroup.id}\\tLB:~{readgroup.lib_id}\\tSM:~{sample.id}\\tPL:~{platform}",
                    bwaIndex = bwaIndex
            }

            # Exclude supplementary alignments
            call samtools.View as filterReads {
                input:
                    inFile = bwaMem.outputBam,
                    outputFileName = readgroupDir + "/" + sample.id + "-" + readgroup.lib_id + "-" + readgroup.id + ".filtered.bam",
                    excludeFilter = 2048
            }
        }
        
        call picard.MarkDuplicates as removeDuplicates {
            input:
                inputBams = filterReads.outputBam,
                outputBamPath = sampleDir + "/~{sample.id}.remove_dups.bam",
                metricsPath = sampleDir + "/~{sample.id}.remove_dups.metrics",
                removeDuplicates = true
        }

        call bammetrics.BamMetrics as metrics {
            input:
                bam = removeDuplicates.outputBam,
                bamIndex = removeDuplicates.outputBamIndex,
                outputDir = sampleDir,
                referenceFasta = referenceFasta,
                referenceFastaFai = referenceFastaFai,
                referenceFastaDict = referenceFastaDict
        }

        call freec.Freec as copyNumberCalling {
            input:
                bamFile = removeDuplicates.outputBam,
                bamIndex = removeDuplicates.outputBamIndex,
                referenceFastaFai = referenceFastaFai,
                mappability = select_first([mappability, gemMappability.mappability]),
                chrFiles = select_first([chrFiles, splitFasta.fastaFiles]),
                outputDir = "~{sampleDir}/FREEC",
                pairedEnd = defined(sample.readgroups[0].R2)
        }

        call freec.makeGraph2_0 as cnvPlots {
            input:
                ratio = copyNumberCalling.ratio,
                outputDir = "~{sampleDir}/FREEC"
        }

        call shallowHRD.ShallowHRD_hg19_controlfreec_chrX as calculateHRD {
            input:
                ratioTxt = copyNumberCalling.ratio,
                outputDir = "~{sampleDir}/shallowHRD"
        }

        Array[File] sampleReports = flatten([metrics.reports,
                                             [removeDuplicates.metricsFile],
                                             flatten(qualityControl.reports)])
        String sampleName = sample.id
    }

    call GetHRDStatus as getHRDStatus {
        input:
            numberLGAs = calculateHRD.number_LGAs,
            samples = sampleName
    }

    Array[File] allReports = flatten(sampleReports)

    call multiqc.MultiQC as multiqcTask {
        input:
            reports = allReports
    }

    output {
        Array[Array[File]] shallowHRD_results = calculateHRD.all
        Array[File] freecGcProfile = copyNumberCalling.gcProfile
        Array[File] freecCnv = copyNumberCalling.cnv
        Array[File] freecInfo = copyNumberCalling.info
        Array[File] freecRatio = copyNumberCalling.ratio
        Array[File] freecSampleCpn = copyNumberCalling.sampleCpn
        Array[File] freecCnvPlots = cnvPlots.ratioPng
        File hrdStatus = getHRDStatus.hrdStatus
        Array[File] filteredBams = removeDuplicates.outputBam
        Array[File] filteredBamIndexes = removeDuplicates.outputBamIndex
        Array[File] reports = allReports
        File multiqcReport = multiqcTask.multiqcReport
    }
}


task SplitFasta {
    input {
        File fasta
        String outDir
    }

    command <<<
        set -e 
        mkdir -p ~{outDir}

        python <<EOF
        with open("~{fasta}", "r") as in_fasta:
            line = next(in_fasta, False)
            while line:
                chr, *_ = line[1:].split()
                out_path = f"~{outDir}/{chr}"
                print(out_path)
                with open(out_path, "w") as out_fasta:
                    out_fasta.write(line)
                    line = next(in_fasta, False)
                    while line and line[0] != ">":
                        out_fasta.write(line)
                        line = next(in_fasta, False)
        EOF
    >>>

    output {
        Array[File] fastaFiles = read_lines(stdout())
    }

    runtime {
        cpu: 1
        memory: "1GiB"
        time_minutes: 120 # !UnknownRuntimeKey
        docker: "python:3.12-slim"
    }

    parameter_meta {
        fasta: {description: "The fasta to split up.", category: "required"}
        outDir: {description: "The directory to write the output to.", category: "required"}
    }
}

task GetHRDStatus {
    input {
        Array[File]+ numberLGAs
        Array[String]+ samples
    }

    command <<<
        python <<EOF > HRD_status.tsv
        paths = ["~{sep='", "' numberLGAs}"]
        samples = ["~{sep='", "' samples}"]

        print("sample\tHRD\tNumber LGAs 10Mb")
        for sample, path in zip(samples, paths):
            with open(path, "r") as num_lga:
                for line in num_lga:
                    size, num = line.strip().split()
                    if size != "10":
                        continue
                    num = int(num)
                    if num < 15:
                        hrd = "No (< 15)"
                    elif num >= 20:
                        hrd = "Yes (>= 20)"
                    else:
                        hrd = "Borderline [15;19]"
                    print(f"{sample}\t{hrd}\t{num}")
                    break
        EOF
    >>>

    output {
        File hrdStatus = "HRD_status.tsv"
    }

    runtime {
        cpu: 1
        memory: "1GiB"
        time_minutes: length(samples) # !UnknownRuntimeKey
        docker: "python:3.12-slim"
    }

    parameter_meta {
        numberLGAs: {description: "The number_LGAs.txt files from shallowHRD.", category: "required"}
        samples: {description: "The sample names.", category: "required"}
    }
}