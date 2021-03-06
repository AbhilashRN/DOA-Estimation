# DOA-Estimation
Various Direction Of Arrival Algorithms written in MATLAB 

Array processing plays an important role in many diverse application areas. Most modern radar and sonar systems rely on antenna arrays or hydrophone arrays as an essential component of the system. Direction of arrival Estimation has become increasingly important and a staple in various multiple array systems and an integral component in communication systems of the future (5G). This project seeks to analyse various Direction of Arrival algorithms and produce a comparative analysis of the performance of these algorithms with respect to various parameters including Number of elements, resolution, Number of snapshots etc.

Various DOA algorithms like Bartlett, MVDR, MUSIC and ESPRIT with their pseudo spectrum equations were analysed and simulated in MATLAB for various conditions like minimum resolution, varying SNR and time varying snapshots antenna systems. These simulations were initially run on a uniform linear array system and gradually extrapolated to a planar array. Adaptive Beamforming was also implemented in the form Sample Matrix Inversion (SMI). 

Comparative Analysis of the various algorithms clearly painted a picture with MUSIC showing higher and sharper peaks over Bartlett and MVDR. Resolution of the N-array setup was also found and visualised which has been grouped in a table along with the results. ESPRIT although inaccurate compared to MUSIC was found to be computationally efficient. Other parameters like SNR and coherence of the signal sources were also simulated. SMI Beamformer is also implemented on MATLAB resulting in an optimised version of MVDR beamformer with nulls at the interfering signal.

This is a compilation of various algorithms written and work done in my prefinal undergrad year with [Abhishek Prakash](https://www.linkedin.com/in/abhi-prakash-a86932140/) and [Dr. Shushrutha](https://rvce.edu.in//ec-ksshushrutha) as our advisor.

An overview of the project along with the results can be viewed [here](Direction_Of_Arrival_Deck.pdf). 
 
