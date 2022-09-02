#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process simulate_data {
    publishDir "simulated_data", mode: "copy"

    input:
    path simulator_xml
    each r
    each m
    each s1
    each s2
    each rep

    output:
    path "simulator.r${r}_m${m}_s1${s1}_s2${s2}_rep${rep}.tree", emit: tree
    path "simulator.r${r}_m${m}_s1${s1}_s2${s2}_rep${rep}.traj", emit: traj 
    val r, emit: r
    val m, emit: m
    val s1, emit: s1
    val s2, emit: s2
    val rep, emit: rep


    script:
    """
    java -jar ~/code/beast_and_friends/remaster/out/artifacts/remaster_jar/remaster.jar \
        -D Re1=1.5 \
        -D Re2=1.1 \
        -D bu=1 \
        -D r=$r \
        -D m=$m \
        -D s1=$s1 \
        -D s2=$s2 \
        -seed $rep \
        -overwrite \
        $simulator_xml
    """
}

process analyze {
    input:
    each path(analysis_xml)
    path tree
    val r
    val m
    val s1
    val s2
    val rep

    output:
    path "output.log", emit: log

    script:
    """
    java -jar ~/code/beast_and_friends/bdmm-prime/out/artifacts/bdmm_prime_jar/bdmm-prime.jar \
        -overwrite \
        -D tree="$tree" \
        -D removalProb="$r" \
        -D outputFile="output.log" \
        $analysis_xml
    """
}

process process_traces {
    input:
    path trace
    val r
    val m
    val s1
    val s2
    val rep

    output:
    path "processed.log", emit: processed_log

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    read_tsv("$trace") %>% slice_tail(prop=0.9) %>%
    mutate(r=$r, m_truth=$m, s1_truth=$s1, s2_truth=$s2, rep=$rep) %>%
    pivot_longer(cols=c(R0Values, sampPropValues)) %>%
    write_csv("processed.log")
    """
    
}

process tabulate {
    publishDir "figures", mode: "copy"

    input:
    val table_file_name
    path trace

    output:
    path "$table_file_name"

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    read_csv("$trace") %>%
    group_by(r, m_truth, s1_truth, s2_truth, types, name) %>%
    summarize(median=median(value),
        lower=quantile(value, probs=0.025),
        upper=quantile(value, probs=0.975)) %>%
    write_tsv("$table_file_name")
    """
}


workflow simulate {
    main:
    simulate_data(Channel.fromPath("XMLs/simulator.xml"),
                  Channel.of(1,0), // r
                  Channel.of(0.3,0.8), // m
                  Channel.of(0.1,0.2), // s1
                  Channel.of(0.1,0.2), // s2
                  1..5) // rep
    emit:
    tree = simulate_data.out.tree
    r = simulate_data.out.r
    m = simulate_data.out.m
    s1 = simulate_data.out.s1
    s2 = simulate_data.out.s2
    rep = simulate_data.out.rep
}

workflow one_type {
    take:
    tree
    r
    m
    s1
    s2
    rep

    main:
    analyze(Channel.fromPath("XMLs/analyze_1type.xml"),
            tree, r, m, s1, s2, rep)

    process_traces(analyze.out.log, r, m, s1, s2, rep)

    tabulate("results_1type.txt",
             process_traces.out.processed_log.collectFile(keepHeader: true, skip: 1))   
}

workflow two_type {
    take:
    tree
    r
    m
    s1
    s2
    rep

    main:
    analyze(Channel.fromPath("XMLs/analyze_2type.xml"),
            tree, r, m, s1, s2, rep)

    process_traces(analyze.out.log, r, m, s1, s2, rep)

    tabulate("results_2type.txt",
             process_traces.out.processed_log.collectFile(keepHeader: true, skip: 1))
}


workflow {
    simulate()

    one_type(simulate.out.tree,
             simulate.out.r,
             simulate.out.m,
             simulate.out.s1,
             simulate.out.s2,
             simulate.out.rep)

    two_type(simulate.out.tree,
             simulate.out.r,
             simulate.out.m,
             simulate.out.s1,
             simulate.out.s2,
             simulate.out.rep)
}
