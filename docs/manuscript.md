---
title: manuscript
titlepage: false
author: ""
date: ""
output:
  docx:
    toc: false
    output: manuscript.docx
    bibliography: /Users/alexander.kvist/Zotero/zotero_library.bib
    csl: /Users/alexander.kvist/Zotero/styles/apa.csl
    citeproc: true
    filter:
      - pandoc-acro
      - pandoc-crossref
acronyms:
  PD:
    short: PD
    long: Parkinson's disease
  fNIRS:
    short: fNIRS
    long: functional near-infrared spectroscopy

---

# Introduction

Text

# Material and methods

## Hypotheses

Introduce hypotheses here? H1-H3.

## Participants

People with Parkinson's disease (≥ 60 years, clinical diagnosis ≥ 6 months prior to enrollment) and age- and sex-matched controls were re-invited from the Park-MOVE cohort [@franzen2023], who had previously participated in an fNIRS study on complex walking. Exclusion criteria for the control group were medical conditions affecting gait and balance, or severe hearing or visual impairment. Exclusion criteria for the Parkinson group included other neurological diseases, severe hearing or visual impairment, other medical conditions affecting gait and balance, and speech difficulties such as aphasia.

## Experimental procedure

Data collection took place at the uMOVE core facility, Karolinska University Hospital, Solna, Stockholm. All data was collected during a single experimental session. Collected data consisted of resting-state fNIRS data, cognitive screening data (MoCA), and a clinical test of balance. The Parkinson group also performed the Movement Disorder Society sponsored revision of the Unified Parkinson’s Disease Rating Scale (MDS-UPDRS) to assess disease severity at the time of measurement.

## Data acquisition

The fNIRS system used was a NIRSport2 (NIRx) with 16 sources and 15 detectors, with 8 short-separation detectors. Optodes transmitted light at 760 and 850 nm. Sampling frequency was 7.6 Hz. Because optodes were not sufficient to ensure full-brain coverage, two measurements were taken in succession with one montage for the left hemisphere, and one for the right hemisphere (Figure 1). Each resting-state measurement lasted for 10 minutes.

The resting-state measurements took place in a calm, dimly lit room, with participants seated in an office chair, legs placed on a leg rest to decrease movement. Eyes were closed, and participants were blindfolded as well. Ear plugs were used to shield from ambient noise. Participants were instructed not to focus on any special thought, let the mind wander, and to not fall asleep.

## Data preprocessing

Preprocessing of resting-state fNIRS data was performed with the MATLAB NIRS AnalyzIR toolbox (forked version; see Data availability for details) [@santosa2018]. Raw data was trimmed at the start and end of the measurement with 10 seconds on each side to avoid noise from position adjustments and similar. Raw data was then converted into optical density and into ΔHbO2 and ΔHHb using the modified Beer-Lambert law [@delpy1988]. The differential path-length factor (DPF) was set to depend on age [@scholkmann2013].

Data quality was assessed using MNE-NIRS (v0.7.1) [@gramfort2013; @luke2021]. The scalp coupling index (SCI) was calculated on a per-channel basis and the percentage of bad channels (SCI < 0.7) among participants was plotted on the montage (supplementary). The worst quality channels were discarded from further analysis, 6 in total, leaving 44 long channels remaining. Moreover, because noise in short-channels may transfer noise into the signal of interest during short-channel regression, quality of the short-channel data was further assessed using their power spectrums [@novi2023] and channels were discarded if they lacked a clear heart-rate peak around 1 Hz (details of discarded channels in supplementary material).

## Data analysis

All-to-all channel correlation matrices were calculated from ΔHbO2 values using the connectivity module in the NIRS toolbox, with prewhitening and short-channel regression [@santosa2017; @lanka2022].

To obtain graph metrics, in particular eigenvector centrality, correlation matrices were fed into the BRAPH2.0 (v2.0.0.b2) [@mijalkov2017] pipeline Connectivity Comparison WU (weighted undirected). This pipeline compares two groups of subjects by constructing weighted undirected graphs from input connectivity data, in this case the correlation matrices. Brain atlases for the BRAPH2.0 pipeline were obtained by exporting the projected MNI coordinates of channels in the used montages via AtlasViewer (v2.44.0) [@aasted2015].

Differences in eigenvector centrality values were compared by non-parametric permutation tests with 10000 permutations, considered significant for a two-tailed t-test at p < .05. To account for multiple comparisons in the network results, false discovery rate (FDR) correction was applied, and the network was tested at q < .05.

Finally, to test hypotheses in H3, correlation between eigenvector centrality and disease severity (MDS-UPDRS motor score), levodopa equivalent daily dose (LEDD), disease duration and MoCA score were calculated. Correlations were calculated for those channels showing the largest differences in eigenvector centrality values between the groups.
