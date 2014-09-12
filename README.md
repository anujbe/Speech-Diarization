Speech-Diarization
==================

Winter Internship-2010-CMU &amp; IIIT-Speech Segmentation &amp; Clustering with Distance-Based Segmentation, BIC and Gaussian Distribution

BIC Segmentation:
  Distance-based segmentation can segment the speech into excessive number of speech change points.
  By undergoing a second pass using BIC method we more accurately identify the speech change points

BIC Clustering:
  BIC clusters based on generalized maximum likelihood criteria. 
    Heirarchical Clustering : This is used to locate all the segments corresponding to the same speaker. 
    Clustering can be performed in a bottom-up manner i.e starting with all the initial segments and building a tree of clusters by merging them

References:
  Audio data indexing : use of second-order statistics for speaker-based segmentation Perrine Delacourt and Christian Wellekens 
  Blind Clustering Of Speech Utterances based on speaker and language characteristics. A. Reynolds, E. Singer, B.A. Carlson 
  Improved Speaker Segmentation and Segments Clustering using the Bayesian Information Criterion Alain Tritschler and Ramesh Gopinath
