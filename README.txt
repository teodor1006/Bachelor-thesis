% --------------------------------------------------------
%      Bachelor thesis-README
%      =========================================
%      File:                 README.txt
%      Author(s):            Teodor Slaveykov
% --------------------------------------------------------

This README presents the content used to create this bachelor's thesis.

- The "Literature" folder contains all images and PDF files used in the context of this work.

- The PDF file "Teodor_BA" presents the entire bachelor's thesis in a readable format.

- The PDF file "Extended Abstract" summarizes the main content and results of the work.

The "DQC_Partition" folder contains various Python files.

:: "GA.py" is responsible for the random partitioning of qubits using the Genetic Algorithm into 4 partitions.

:: "KerLin.py" considers qubit partitioning using the Kernighan-Lin algorithm into 4 partitions.

:: "Spectral.py" focuses on Spectral Partitioning into 4 partitions.

:: "GA_MKL.py" and "GA_SP.py" are based on 2 hybrid Genetic Algorithms (HGA).

:: "GA_KL.py" verifies the results based on the scientific article by Zomorodi & Co., which considers only
2 partitions.

:: "QiskitTranspiler.py" converts the benchmark circuits from RevLib into OpenQasm 2.0 format.

:: "Experiments_for_Verification.py" presents and compares the results of Zomorodi & Co. with mine for 2
partitions.

:: The "qasm" folder contains 9 benchmark circuits in OpenQasm 2.0 format.




