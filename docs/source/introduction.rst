

Introduction to spike-in control
--------------------------------
`ChIP-sequencing (ChIP-Seq) <https://en.wikipedia.org/wiki/ChIP_sequencing>`_
is widely used to map transcription factor binding sites (TFBS) and histone modifications (methylation, acetylation, phosphorylation, and ubiquitylation) status to the genome. One important step of ChIP-seq data analysis is **peak calling**, in which the ChIPed sample is compared to the control sample to find peak regions where "DNA fragments from ChIPed samples" are significantly enriched.

Besides **peak calling**, ChIP-seq data can also be used for **quantitative analyses**, in which one group of ChIPed samples are compared to another group of ChIPed samples to identify
differential TFBS, gained or lost histone modification sites, etc.

Both **peak calling** and **quantitative analysis** assume the same number of cells would yield the same amount of DNAs (or chromatin) under different conditions. Sometimes, this assumption does not hold [1]_. For example, when the histone methyltransferases (such as EZH2 for H3K27 [2]_, SETD2 for H3K36 [3]_, DOT1L for H3K79 [4]_) are mutated or dysregulated. The per-cell DNA/chromatin yield in ChIP experiments will be significantly reduced. In these scenarios, conventional normalization and analytic methods will fail or compromise the power to identify peaks or detect epigenetic changes.

To accurately quantify the ChIP-seq signals, we need to set up an "internal control"--a small amount of exogenous chromatin (called spike-in chromatin) added to all the samples before doing the ChIP experiment. Any biological manipulations or technical variations introduced during the ChIP-seq experiment will also occur with the spike-in chromatin. After sequencing, we can normalize the signals from the experimental samples to this "internal control". 

The `Drosophila melanogaster <https://en.wikipedia.org/wiki/Drosophila_melanogaster>`_ (fruit fly) genome is a good exogenous control for mammalian cells because 1) The Drosophila genome is well studied and has a high-quality sequence assembly. 2) The similarity between Drosophila and human (mouse) genomes is low, so that it's easy to separate "human DNA sequences" from "Drosophila DNA sequences". 3) Drosophila cells are readily available in large quantities. 4) The Drosophila cells have all of the key histone modification marks reported in humans. 


**References**


.. [1] Chen K, Hu Z, Xia Z, Zhao D, Li W, Tyler JK. The Overlooked Fact: Fundamental Need for Spike-In Control for Virtually All Genome-Wide Analyses. Mol Cell Biol. 2015 Dec 28;36(5):662-7. doi: 10.1128/MCB.00970-14. PMID: `26711261 <https://pubmed.ncbi.nlm.nih.gov/26711261/>`_

.. [2] Cao R, Wang L, Wang H, Xia L, Erdjument-Bromage H, Tempst P, Jones RS, Zhang Y. Role of histone H3 lysine 27 methylation in Polycomb-group silencing. Science. 2002 Nov 1;298(5595):1039-43. doi: 10.1126/science.1076997. Epub 2002 Sep 26. PMID: `12351676 <https://pubmed.ncbi.nlm.nih.gov/12351676/>`_.

.. [3] Pfister SX, Ahrabi S, Zalmas LP, Sarkar S, Aymard F, Bachrati CZ, Helleday T, Legube G, La Thangue NB, Porter AC, Humphrey TC. SETD2-dependent histone H3K36 trimethylation is required for homologous recombination repair and genome stability. Cell Rep. 2014 Jun 26;7(6):2006-18. doi: 10.1016/j.celrep.2014.05.026. Epub 2014 Jun 12. PMID: `24931610 <https://pubmed.ncbi.nlm.nih.gov/24931610>`_

.. [4] Orlando DA, Chen MW, Brown VE, Solanki S, Choi YJ, Olson ER, Fritz CC, Bradner JE, Guenther MG. Quantitative ChIP-Seq normalization reveals global modulation of the epigenome. Cell Rep. 2014 Nov 6;9(3):1163-70. doi: 10.1016/j.celrep.2014.10.018. Epub 2014 Oct 30. PMID: `25437568 <https://pubmed.ncbi.nlm.nih.gov/25437568/>`_.
