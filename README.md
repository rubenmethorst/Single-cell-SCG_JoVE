# Single-cell-SCG_JoVE
Authors: Y. Ge, L. van Roon, H. Chen, R. Methorst, M. Paton, M.C. de Ruiter, S.M. Kielbasa, M.R.M. Jongbloed

This repository holds all scripts used to analyze data and create figures presented in: Low-input nucleus isolation and multiplexing with barcoded antibodies of mouse sympathetic ganglia for single-nucleus RNA sequencing 

# **Abstract**
The cardiac autonomic nervous system is crucial in controlling cardiac function, such as heart rate and cardiac contractility, and is divided into a sympathetic and parasympathetic branch. In healthy persons there is a balance between these two branches to maintain homeostasis. However, cardiac disease states such as myocardial infarction, heart failure and hypertension can induce remodeling of cells involved in cardiac innervation, which in turn is associated with adverse clinical outcome. Although a myriad of data is present on the histological structure and function of the cardiac autonomic nervous system, its molecular biological architecture in health and disease is in many aspects still enigmatic. Novel technologies such as single cell RNA sequencing hold promise for genetic characterization of tissues on a single cell level, however the relatively large size of neurons may impede the standardized use of these techniques. Here, This protocol exploits droplet-based single-nucleus RNA sequencing (snRNA-seq), a promising method to characterize the biological architecture of cardiac sympathetic neurons in health and in disease. A stepwise approach is demonstrated to perform snRNA-seq of the bilateral superior cervical and stellate ganglia dissected from adult mice. This method enables long-term sample preservation maintaining an adequate RNA quality when samples cannot be fully collected within a short period of time. Nucleus-barcoding with hashtag barcoding antibody-oligos (HTOs) staining enables demultiplexing and the trace-back of distinct ganglionic samples during the afterward single nucleus analysis. The analysis results support successful nuclei capture of neuronal cells, glial cells and endothelial cells of the sympathetic ganglion by snRNA-seq. In summary, this protocol provides a stepwise approach for single nucleus RNA sequencing of sympathetic cardiac ganglia, a method that has the potential for broad application in studies of innervation of other organs and tissues.

# **Script**

* `Script_SCG_JOVE.R`: Used to perform all analyses presented in the manuscript

Figures can also be plotting using above mentioned script. Associated results are located here: `Result/`
In the `Data/` folder, there is a excel-file containing a variety of cell markers for visualization purposes


--------------------------
## **The MIT License (MIT)**

### Copyright (c) 1979-2020

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.

