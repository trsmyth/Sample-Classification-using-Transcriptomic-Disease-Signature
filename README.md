# Sample-Classification-using-Transcriptomic-Disease-Signature

This repository contains code for a personal side project seeking to classify human tissue samples according to transcriptomic signature based on publicly available disease signature designations.

This code seeks to identify primary human tissue samples from the Disease Signature Atlas (DiSignAtlas), assign appropriate metadata, access uniformly processed and aligned count data, and assign each sample a corresponding disease signature based on Medical Subject Heading (MeSH) designations. 

Cleaned metadata and count data is currently being used to build a deep neural net (DNN) using Pytorch with the aim of developing a classifier for each control tissue type and MeSH desingnation (ongoing development).

DiSignAtlas : http://www.inbirg.com/disignatlas/download (Gene Expression Profiles of Each DiSignAtlas Dataset)
ARCHS4 version 2.5 : https://maayanlab.cloud/archs4/download.html
MeSH designations : https://nlmpubs.nlm.nih.gov/projects/mesh/MESH_FILES/mesh_browser/ (desc2024, supp2024)
