**JAYPEE INSTITUTE OF INFORMATION TECHNOLOGY**

**B.Tech Major Project – Semester VII**

**PROJECT SYNOPSIS**

**Performance Analysis of OTFS under High-Mobility Multipath and Fading Channels**

**Batch: A17**

| Team Member 1 | Arnav Maindola – 23119021 |
| :---- | :---- |
| Team Member 2 | Mukta Pande – 23119019 |
| Supervisor | Prof. Alok Joshi |
| Programme | B.Tech Electronics & Communication Engineering |

# **1\. Introduction**

High-mobility wireless communication environments are characterized by time-varying multipath channels and Doppler shifts, which can degrade the reliability of conventional multicarrier communication systems. Orthogonal Frequency Division Multiplexing (OFDM), although widely adopted because of its spectral efficiency and simple implementation, can suffer from loss of subcarrier orthogonality and increased inter-carrier interference under significant Doppler conditions.

Orthogonal Time Frequency Space (OTFS) modulation has emerged as a promising waveform for high-mobility communication because it represents information in the delay-Doppler domain. This project investigates the performance of OTFS under realistic multipath and different statistical fading conditions, while also comparing OTFS with conventional OFDM under the same high-mobility channel conditions.

# **2\. Problem Statement**

To investigate the robustness of OTFS in high-mobility wireless communication by comparing its performance with OFDM over an Extended Vehicular A (EVA) multipath channel and by evaluating OTFS under different fading environments using BER as the primary performance metric.

# **3\. Aim**

To develop a MATLAB-based simulation framework for analyzing OTFS performance under high-mobility multipath channels and different fading conditions, and to establish its relative robustness compared with OFDM.

# **4\. Objectives**

* Implement an OTFS communication system in MATLAB using QPSK modulation and MMSE-based detection.  
* Model a realistic Extended Vehicular A (EVA) multipath channel with additive white Gaussian noise (AWGN).  
* Compare OFDM and OTFS under identical EVA channel and mobility conditions using BER versus SNR and BER versus velocity.  
* Investigate OTFS performance under Rayleigh, Rician, Nakagami-m, α-μ, and κ-μ fading models.  
* Study the effect of increasing mobility and the corresponding Doppler shift on OTFS reliability.  
* Analyze and compare BER trends to identify the channel and mobility conditions under which OTFS provides improved robustness.

# **5\. Proposed Methodology**

The project is divided into two major parts.

## **Part A – OFDM versus OTFS**

OFDM and OTFS systems will be implemented using the same QPSK modulation, EVA multipath channel, AWGN level, SNR range, velocity, and channel realization wherever applicable. The purpose is to provide a controlled comparison in which the waveform is the primary variable. The performance will be evaluated using:

* BER versus SNR  
* BER versus velocity (with velocity producing the corresponding Doppler shift)

## **Part B – OTFS under Different Fading Conditions**

The OTFS system will be evaluated under five statistical fading models: Rayleigh, Rician, Nakagami-m, α-μ, and κ-μ. The fading analysis will be performed in conjunction with the multipath/noise model, and the resulting system performance will be evaluated using:

* BER versus SNR  
* BER versus velocity / Doppler-induced mobility

# **6\. Simulation Specifications**

| Parameter | Specification |
| :---- | :---- |
| Simulation platform | MATLAB |
| Waveforms | OFDM and OTFS |
| Modulation | QPSK |
| OTFS grid | M \= 32, N \= 16 |
| Subcarrier spacing | 15 kHz |
| Carrier frequency | 4 GHz |
| Primary detector/equalizer | MMSE |
| SNR range | 0–30 dB, 5 dB step |
| Velocity range | 0–500 km/h |
| Multipath profile | Extended Vehicular A (EVA) |
| Noise model | AWGN |
| Fading models in Part B | Rayleigh, Rician, Nakagami-m, α-μ, κ-μ |
| Monte Carlo simulation | Approximately 500–1000 frames for final results |

# **7\. Performance Metrics and Expected Plots**

BER will be the primary performance metric. Four core result plots are planned:

* Part A: BER versus SNR – OFDM versus OTFS over EVA \+ AWGN.  
* Part A: BER versus velocity – OFDM versus OTFS over EVA \+ AWGN.  
* Part B: BER versus SNR – OTFS under the five selected fading models.  
* Part B: BER versus velocity – OTFS under the five selected fading models.

# **8\. Expected Outcomes**

* A validated MATLAB simulation framework for OTFS under high-mobility multipath conditions.  
* A controlled comparison of OFDM and OTFS under an EVA vehicular channel.  
* Quantitative understanding of the effect of SNR and mobility/Doppler on OTFS BER.  
* A comparative assessment of OTFS under Rayleigh, Rician, Nakagami-m, α-μ, and κ-μ fading.  
* Identification of channel conditions in which OTFS demonstrates improved robustness compared with OFDM.

# **9\. Significance of the Work**

The study provides a systematic evaluation of OTFS in conditions relevant to high-mobility wireless communication. The combination of a standardized vehicular multipath profile and multiple statistical fading models enables the limitations and robustness of OTFS to be studied from both realistic and generalized channel perspectives. The simulation framework can also serve as the foundation for further waveform and propagation-enhancement studies.

# **10\. Future Scope – Semester VIII**

The validated OTFS framework developed in Semester VII can be extended in Semester VIII towards intelligent reflecting surface (IRS)-assisted OTFS communication. The extension will investigate whether controlled reflected propagation can improve OTFS performance under the challenging fading and high-mobility conditions identified in this study.

# **11\. Conclusion**

The proposed Semester VII project focuses on establishing the high-mobility performance of OTFS through a controlled OFDM comparison and a detailed fading analysis. By using EVA multipath propagation, AWGN, multiple fading distributions, QPSK modulation, and MMSE detection, the study aims to provide a reproducible and technically meaningful assessment of OTFS for high-mobility wireless communication.