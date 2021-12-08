# MFT_CPU_Tutorial
Simple Tutorial for MFT (CPU Version), source codes courtesy Dr. Meng (SCEC)


a. Continuous waveform 
Sorted under a directory ContWF/, named as each day: 20190610, for instance. 
Put all waveforms for used stations together: make sure all of them have roughly
the same length (the problem would use common window for all traces); all of them
shift back to the same origin time (no negative values allowed). 

About the filter: this generally depends on the purpose. To filter out signals from 
regional/teleseismic events as well as suppress culture (traffic) noises, a proper 
filter would be 2–8 Hz or 2–16 Hz. If no high frequency content is needed, it would 
save a lot of computation by downsampling the data (to 20 for 2–8 Hz or 40 for 2–16 Hz)

b. Template waveform 
To avoid potential mismatching by pre-processing data (cutting from raw data, filtering, 
and downsampling), it would cause a non-one cross-correlation between the "same" segment 
of data from requesting individually or cutting from continuous trace.

After cutting waveform from continous data, phase picking could be time-consuming. A 
good candidate would be a certain kind of phase picker (obspy picker, for instance). 
Then calculate the SNR for all phase picks (to determine whether this phase window 
could be used as a matched filter). 

c. Working direcotry 
