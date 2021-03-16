# required modules R subread

import glob, os, yaml, pandas

maindir = os.path.abspath(workflow.basedir)

print(config)

samplesheet_files, = glob_wildcards(os.path.join(maindir,config['samplesheets_dir'],'{samplesheet}.tsv'))

# parameter: log2fold cutoff and padj in APA target identification
assert all(x in config.keys() for x in ['padj_cutoff']) , 'Please make sure that you define a padj cutoff for APA target genes'

samplesheets = dict()
for n in samplesheet_files:
    samplesheets[n] = os.path.join(maindir, config['samplesheets_dir'],''.join([n,'.tsv']))

bam_files = glob.glob(os.path.join(maindir,config['bam_dir'],'*.bam'))
print(bam_files)

paqr_files = ','.join(config['paqr'])
isoscm_files = ','.join(config['isoscm'])

rule all:
    input:
        os.sep.join([config['project_name'],"config.yaml"]),
        expand(os.sep.join([config['project_name'],'APA_targets',"{samplesheet}.APA_targets.tsv"]), samplesheet = samplesheets.keys())


rule writeConfig:
    output:
        yaml=os.sep.join([config['project_name'],"config.yaml"])
    run:
        yaml_stream = open(output.yaml,'w')
        yaml.dump(config, yaml_stream)
        yaml_stream.close()

# potentially outsource and provide as external tool, if the user has isoSCM and PAQR
rule poolBreakpoints:
    input:
        annotation=config['annotation'],
        segments=config['segments']
    output:
        breakpoints=os.sep.join([config['project_name'], 'Annotation',"Breakpoints_pooled.merged_downstream_breakpoint.gff"])
    params:
        paqr=','.join(config['paqr'].split(' ')),
        isoscm=','.join(config['isoscm'].split(' ')),
        min_dist=config['min_distance'],
        isoscm_confidence=config['isoscm_confidence'],
        prefix=os.sep.join([config['project_name'], 'Annotation',"Breakpoints_pooled"])
    log:
        os.path.join(config['project_name'], 'Annotation','log', 'Breakpoints_pooled.log')
    shell:
        "Rscript {script} {options} {gtf} {segments} {isoscm_files} {paqr_files} {outfile} &> {log}".format(
            script = os.sep.join([maindir, "tools/mergeBreakpoints.R"]),
            options = "--isoscm.confidence={params.isoscm_confidence} --minimum.distance={params.min_dist}",
            gtf="--gtf {input.annotation}",
            segments="--segments {input.segments}",
            isoscm_files="--isoscm {params.isoscm}",
            paqr_files="--paqr {params.paqr}",
            outfile="--outfileNamePrefix {params.prefix}",
            log = "{log}")

# Generalize to use any polyA database
# make optional
rule breakSegments:
    input:
        breakpoints=rules.poolBreakpoints.output.breakpoints,
        segments=config['segments'],
        annotation=config['annotation']
    output:
        nodes_split=os.sep.join([config['project_name'], 'Annotation',"Segments_split.gff"]), # <prefix>.gff
        saf=os.sep.join([config['project_name'], 'Annotation',"Segments_split.saf"])          # <prefix>.saf
    params:
        min_dist=config['min_distance'],
        isoscm_confidence=config['isoscm_confidence']
    log:
        os.path.join(config['project_name'], 'Annotation', 'log', 'Segments_split.log')
    shell:
        "Rscript {script} {options} {gtf} {segments} {breakpoints} {outfile} --debug &> {log}".format(
            script = os.sep.join([maindir, "tools/breakWhippetNodes.R"]),
            options = "--isoscm.confidence={params.isoscm_confidence} --min.distance={params.min_dist}",
            gtf="--gtf {input.annotation}",
            segments="--segments {input.segments}",
            breakpoints="--breakpoints {input.breakpoints}",
            outfile="--outfileName {output.nodes_split}",
            log = "{log}")

rule featureCounts:
    input:
        bams=bam_files,
        saf=rules.breakSegments.output.saf
    output:
        tsv=os.sep.join([config['project_name'], 'featureCount',''.join([config['project_name'],".segments_split.featureCounts.tsv"])])
    params:
        # -f per feature
        # -s 2 library type
        # --primary reads only
        # -p paired
        # -B both aligned
        others="-t sequence_feature -g feature_id -f -s 2 --primary -p -B"
    threads: 16
    log:
        os.path.join(config['project_name'], 'featureCount','log',''.join([config['project_name'],".segments_split.featureCounts.log"]))
    shell:
        """
        featureCounts {params.others} -T {threads} \
            -a {input.saf} -F SAF -o {output.tsv} {input.bams} &> {log}
        """

rule Quantify:
    input:
        featureCounts=rules.featureCounts.output.tsv,
        samplesheet=lambda wildcard: samplesheets[wildcard.samplesheet]
    output:
        dexseq=os.sep.join([config['project_name'], 'DEXseq',"{samplesheet}.segments_split.dexseq.tsv"])
    log:
        os.sep.join([config['project_name'], 'DEXseq','log',"{samplesheet}.segments_split.dexseq.log"])
    shell:
        "Rscript {script} {featurecounts} {samplesheet} {output} &> {log}".format(
            script = os.sep.join([maindir, "tools/DEXseq_fromWhippetFeatureCounts.R"]),
            featurecounts="--table {input.featureCounts}",
            samplesheet="--samplesheet {input.samplesheet}",
            output="--outfileName {output.dexseq}",
            log = "{log}")

rule APA_targets:
    input:
        annotation=config['annotation'],
        segments=rules.breakSegments.output.nodes_split,
        dexseq=rules.Quantify.output.dexseq
    output:
        tsv=os.sep.join([config['project_name'], 'APA_targets',"{samplesheet}.APA_targets.tsv"]),
    params:
        padj_cutoff=config['padj_cutoff'],
        dexseq2gff='--DEXseq2gff',
        apa2gff='--APA2gff',
        write_locus='--write_locus',
        width_min ="--width.min 100"
    log:
        os.path.join(config['project_name'], 'APA_targets','log',"{samplesheet}.APA_targets.log")
    shell:
        "Rscript {script} {annotation} {segments} {dexseq} {width_min} {options} {output} &> {log}".format(
            script = os.sep.join([maindir, "tools/APA_target_identification.R"]),
            annotation="--gtf {input.annotation}",
            segments="--segments {input.segments}",
            dexseq="--dexseq.table {input.dexseq}",
            width_min="{params.width_min}",
            options=" --padj.cutoff {params.padj_cutoff} {params.dexseq2gff} {params.apa2gff} {params.write_locus}",
            output="--outfileName {output.tsv}",
            log = "{log}")
