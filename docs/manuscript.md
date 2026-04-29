---
title: manuscript
titlepage: false
author: ""
date: ""
output:
  docx:
    toc: false
    output: manuscript.docx
    metadata:
      csl: /Users/alexander.kvist/Zotero/styles/apa.csl
    bibliography: /Users/alexander.kvist/Zotero/zotero_library.bib
    filter:
      - pandoc-acro
      - pandoc-crossref
    citeproc: true
    verbose: true
acronyms:
  PD:
    short: PD
    long: Parkinson's disease
  AD:
    short: AD
    long: Alzheimer's disease
  fNIRS:
    short: fNIRS
    long: functional near-infrared spectroscopy
  fMRI:
    short: fMRI
    long: functional magnetic resonance imaging
  BOLD:
    short: BOLD
    long: blood-oxygen level dependent response
  rsFC:
    short: rsFC
    long: resting-state functional connectivity
  IPL:
    short: IPL
    long: inferior parietal lobule
  STG:
    short: STG
    long: superior temporal gyrus
  SMG:
    short: SMG
    long: supramarginal gyrus
  PFC:
    short: PFC
    long: prefrontal cortex
  HbO:
    short: HbO
    long: oxygenated hemoglobin
  HbR:
    short: HbR
    long: deoxygenated hemoglobin
  MoCA:
    short: MoCA
    long: Montreal Cognitive Assessment
  LEDD:
    short: LEDD
    long: levodopa equivalent daily dose
  OSF:
    short: OSF
    long: Open Science Framework

---

# Introduction

When the brain is at rest, it has been shown that there are synchronized low-frequency oscillations of neuronal activity, which can be measured via the +BOLD [@biswal1995]. These oscillations exhibit consistent structural patterns across individuals and are referred to as +rsFC [@damoiseaux2006]. Several +rsFC networks are known in humans and are involved in a wide variety of functions, from memory and executive functioning to motor function [@damoiseaux2006; @uddin2019], and have been extensively studied using +fMRI [@biswal2025].

For some neurodegenerative diseases, these metrics have shown some promise as a biomarker: for +AD, alterations in the default mode network have been found, and in +PD, alterations in limbic and motor related networks have been found [@hohenfeld2018]. For +PD, a meta-analysis of a range of resting-state +fMRI metrics found a convergence of abnormal connectivity patterns converging in the bilateral posterior +IPL, with the direction of effect compared to healthy controls depending on dopaminergic medication [@tahmasian2017]. Specifically, metrics were found to increase compared to controls while off medication, and decrease while on medication [@tahmasian2017], something also found in a study comparing eigenvector centrality between +PD and controls, where eigenvector centrality was found to be lower in the +STG, the +SMG, and the putamen when on medication [@ballarini2018].

While +fMRI remains the most established method for studying resting-state networks, there has been increased investigation into the possibility of studying +rsFC using +fNIRS [@albrecht2025; @wang2024; @yeung2020], an optical neuroimaging method based on detecting changes on +HbO and +HbR concentrations in the cortex [@jobsis1977]. Progress has been made in characterizing the test-retest reproducibility of +rsFC metrics measured with +fNIRS [@zhang2011; @zhang2011a], identifying optimal processing methods [@novi2023; @lanka2022; @paranawithana2022; @blanco2018; @santosa2017], and characterizing resting-state networks [@abdalmalak2022; @butters2026]. Some studies have compared +rsFC metrics and networks identified with +fNIRS to those measured with +fMRI [@kotsogiannis2026; @uchitel2022; @duan2012], with two studies finding a large degree of overlap between network topological parameters and network maps [@uchitel2022; @duan2012], and another finding less agreement for finer-scale nodal parameters but better agreement for large-scale parameters such as global connectivity patterns and major functional subnetworks [@kotsogiannis2026].

However, these comparisons with +fMRI literature have focused on healthy participants, and because resting-state +fNIRS studies are still very heterogeneous  [@albrecht2025], there is a lack of rigorous comparison between resting-state alterations reported in the literature for neurodegenerative disease, and +rsFC metrics in measured with +fNIRS. Thus, this pre-registered study set out to test a number of hypotheses based on the +fMRI literature on resting-state alterations in +PD, in a cohort of +PD participants measured with resting-state +fNIRS.

# Material and methods

## Participants

People with Parkinson's disease (≥ 60 years, clinical diagnosis ≥ 6 months prior to enrollment) and age- and sex-matched controls were re-invited from the Park-MOVE cohort [@franzen2023], who had previously participated in an fNIRS study on complex walking. Exclusion criteria for the control group were medical conditions affecting gait and balance, or severe hearing or visual impairment. Exclusion criteria for the Parkinson group included other neurological diseases, severe hearing or visual impairment, other medical conditions affecting gait and balance, and speech difficulties such as aphasia. Participants in the Parkinson group were on their usual medication schedule.

## Pre-registration and hypotheses

A pre-registration of hypotheses based on +fMRI literature on resting-state alterations in +PD was made at the start of data collection, and is available at +OSF (<https://osf.io/wtsy7/overview?view_only=a8b2c4e92aa440cc87c9448d2e03bd04>). No deviations from the pre-registration are noted.

Based on the literature, we hypothesized that there is a difference in eigenvector centrality between +PD and controls [@ballarini2018] (H1). Furthermore, we expected a lower eigenvector centrality in the +PD group compared to controls in channels most closely corresponding to the MNI coordinates (48 17 -13) in the +STG and (-63 -25 38) in the +SMG [@ballarini2018] (H2.1), as well as MNI coordinates (46 -64 26) in the +IPL [@tahmasian2017] (H2.2).

Finally, we hypothesized that eigenvector centrality in the +PD group is negatively correlated with disease severity, +LEDD, disease duration, and +MoCA score (H3).

## Experimental procedure

Data collection took place at the uMOVE core facility, Karolinska University Hospital, Solna, Stockholm. All data was collected during a single experimental session. Collected data consisted of resting-state fNIRS data, cognitive screening data (MoCA), and a clinical test of balance. The Parkinson group also performed the Movement Disorder Society sponsored revision of the Unified Parkinson’s Disease Rating Scale (MDS-UPDRS) to assess disease severity at the time of measurement.

## Data acquisition

The fNIRS system used was a NIRSport2 (NIRx) with 16 sources and 15 detectors, with 8 short-separation detectors. Optodes transmitted light at 760 and 850 nm. Sampling frequency was 7.6 Hz. Because optodes were not sufficient to ensure full-brain coverage, two measurements were taken in succession with one montage for the left hemisphere, and one for the right hemisphere ([Figure @fig:rs-montage]). Each resting-state measurement lasted for 10 minutes.

The resting-state measurements took place in a calm, dimly lit room, with participants seated in an office chair, legs placed on a leg rest to decrease movement. Eyes were closed, and participants were blindfolded as well. Ear plugs were used to shield from ambient noise. Participants were instructed not to focus on any special thought, let the mind wander, and to not fall asleep.

![Montage used for the resting-state fNIRS measurement, showing optode configuration and sensitivity profiles generated in AtlasViewer for left and right hemispheres](../results_figures/rs_montage.png){#fig:rs-montage}

## Data preprocessing

Preprocessing of resting-state fNIRS data was performed with the MATLAB NIRS AnalyzIR toolbox (forked version; see Data availability for details) [@santosa2018]. Raw data was trimmed at the start and end of the measurement with 10 seconds on each side to avoid noise from position adjustments and similar. Raw data was then converted into optical density and into ΔHbO2 and ΔHHb using the modified Beer-Lambert law [@delpy1988]. The differential path-length factor (DPF) was set to depend on age [@scholkmann2013].

Data quality was assessed using MNE-NIRS (v0.7.1) [@gramfort2013; @luke2021]. The scalp coupling index (SCI) was calculated on a per-channel basis and the percentage of bad channels (SCI < 0.7) among participants was plotted on the montage (supplementary). The worst quality channels were discarded from further analysis, 6 in total, leaving 42 long channels remaining. Moreover, because noise in short-channels may transfer noise into the signal of interest during short-channel regression, quality of the short-channel data was further assessed using their power spectrums [@novi2023] and channels were discarded if they lacked a clear heart-rate peak around 1 Hz ([Figure @fig:short-ch-example]) (details of discarded channels in supplementary material).

![Example of filtering via power spectrum inspection, where the S9-D20 short channel lacks a clear peak around 1 Hz and is therefore discarded.](../results_figures/example_short_ch_check.png){#fig:short-ch-example}

## Data analysis

All-to-all channel correlation matrices were calculated from ΔHbO2 values using the connectivity module in the NIRS toolbox, with prewhitening and short-channel regression [@santosa2017; @lanka2022].

To obtain graph metrics, in particular eigenvector centrality, correlation matrices were fed into the BRAPH2.0 (v2.0.0.b2) [@mijalkov2017; @chang2025] pipeline Connectivity Comparison WU (weighted undirected). This pipeline compares two groups of subjects by constructing weighted undirected graphs from input connectivity data, in this case the correlation matrices. Brain atlases for the BRAPH2.0 pipeline were obtained by exporting the projected MNI coordinates of channels in the used montages via AtlasViewer (v2.44.0) [@aasted2015].

Differences in eigenvector centrality values were compared by non-parametric permutation tests with 10000 permutations, considered significant for a two-tailed t-test at p < .05. To account for multiple comparisons in the network results, false discovery rate (FDR) correction was applied, and the network was tested at q < .05.

To test correlation hypotheses in H3, correlation between eigenvector centrality and disease severity (MDS-UPDRS motor score), +LEDD, disease duration and +MoCA score were calculated. Correlations were calculated for those channels showing the largest differences in eigenvector centrality values between the groups.

Finally, as an exploratory analysis in addition to the eigenvector centrality comparison, the analysis pipeline was run for all other available global and nodal BRAPH graph metrics and compared between groups. 

# Results

## Participants

In total, 57 participants took part in the study (Table 1).

| Variable               | Control (N=30)     | Parkinson (N=27)   | p value |
|------------------------|---------------------|----------------------|---------|
| Gender, female         | 10 (34.5%)          | 10 (37.0%)           | 0.84    |
| Age, yrs               | 68.03 (6.74)        | 69.63 (6.77)         | 0.53    |
| Weight, kg             | 75.31 (13.98)       | 76.81 (19.05)        | 0.99    |
| Height, cm             | 175.24 (10.25)      | 174.63 (9.22)        | 0.79    |
| Mini-BESTest score     | 24.47 (2.43)        | 19.41 (5.09)         | < 0.01  |
| MDS-UPDRS III motor score          | NA                  | 31.78 (14.89)        |        |
| LEDD                   | NA                  | 714.07 (325.12)      |        |
: Table 1 - Demographics

## Eigenvector centrality

Differences in eigenvector centrality values between the groups ([Figure @fig:eigenvector-diff]a, supplementary table 1 and 2) revealed a single channel with a significantly lower value in the Parkinson group compared to the control group (SRC10-DET9, HC 0.153, PD 0.119, difference −0.034, p .034). The center of this channel was according to AtlasViewer simulations located at MNI coordinates (-39 -54 46) ([Figure @fig:eigenvector-diff]b). This location corresponded to Area_PFm_Complex_L in the HCPEx atlas, which is located in the +IPL. However, compared to hypothesized difference in the right +IPL (H2.2), it was located in the left +IPL ([Figure @fig:location-comparison]). Furthermore, the channel did not survive FDR-adjustment.

There were no identified correlations between eigenvector centrality and any of the hypothesized correlates disease severity, +LEDD, disease duration, or +MoCA ([Figure @fig:correlation]).

![(A) Differences in eigenvector centrality values plotted on BRAPH template human_ICBM152. Size and color of circles indicate magnitude and direction of the differences. Blue indicates Parkinson < Control. Red indicates Parkinson > Control. Significant difference from permutation test (p < .05) marked in green. (B) Location of largest difference of eigenvector centrality between the Parkinson group and the control group, visualized on the HCPEx atlas.](../results_figures/braph_diff_combined.png){#fig:eigenvector-diff}

![Comparison between location of largest difference in eigenvector centrality from permutation test (p < .05) (S10-D9) and hypothesized locations of differences in the supramarginal gyrus (SMG) and superior temporal gyrus (STG) (H2.1) and inferior parietal lobule (IPL) (H2.2).](../results_figures/h2_location_comparison.png){#fig:location-comparison}

![Correlation between eigenvector centrality and disease severity (MDS-UPDRS motor score), levodopa equivalent daily dose (LEDD), disease duration and MoCA score for channel SRC10-DET9 in the Parkinson group.](../results_figures/correlation_fig.png){#fig:correlation}

## Exploratory analysis

Exploratory analyses of other graph metrics revealed one channel in the left hemisphere ([Figure @fig:core-periphery-diff]a) that remained after FDR-correction with a higher core-periphery nodal value (SRC3-DET1, difference 0.19, FDR p < .02). The channel was located at MNI coordinates (-21, 52, 1) ([Figure @fig:core-periphery-diff]b). This location corresponded to Area_posterior_10p_L in the HCPEx atlas, located in the left +PFC.

![(A) Differences in core-periphery values plotted on BRAPH template human_ICBM152. Size and color of circles indicate magnitude and direction of the differences. Blue indicates Parkinson < Control. Red indicates Parkinson > Control. Significant difference from permutation test (p < .05) marked in green. (B) Location of largest difference of core-periphery between the Parkinson group and the control group, visualized on the HCPEx atlas.](../results_figures/braph_core_periphery_combined.png){#fig:core-periphery-diff}

# Discussion



# Data availability

All code is available via <https://github.com/alkvi/fnirs_resting_state_pd> and the open lab notebook at <https://alkvi.github.io/fnirs_resting_state_pd>. The original data are not publicly available due to Swedish/EU law, but is located with restricted access in a central repository (<https://doi.org/10.48723/vscr-eq07>), where sharing will be regulated via a data transfer and user agreement upon a reasonable request.

# References
