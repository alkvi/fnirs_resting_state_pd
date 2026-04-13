# Resting-state fNIRS in Parkinson's disease

....

Steps

- Run dataset_preparation/mne_prepare_bids and copy_empty_trigger
- Run quality_control
- Run subject_correlation for left and right
- Copy over left/right R matrices to braph_data HC_left and PD_left etc

- Run AtlasViewer simulation and get MNI coords to create BRAPH atlas, excluding bad channels and shorts

- Install [pandoc-acro](https://github.com/kprussing/pandoc-acro/tree/trunk) in docs/venv

Build manuscript with for example

~~~
pandoc manuscript.md --output manuscript.docx --metadata csl=/Users/alexander.kvist/Zotero/styles/apa.csl --citeproc --bibliography /Users/alexander.kvist/Zotero/zotero_library.bib --filter pandoc-acro --verbose
~~~

or with panrun

~~~
panrun manuscript.md
~~~
