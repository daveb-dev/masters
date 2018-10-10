Rats
For each rat I have saved the following as a matlab variable
(1)time  (in days)
(2)ADC   (apparent diffusion coefficient 4D matrix, 4thdimension is time point)
(3)cells (estimate of the number of tumor cells in 4D, 4thdimension is time point)
(4)brain (mask of everything that is within the brain)
(5)skull (outline of just the skull
(6)anatomical (anatomical often T2w image)
(7)anatomical_full (non-cropped version of the anatomical image)
(8)tissue_segementation
	a.Segmentation of different regions within the brain 
	b.1 = cerebral cortex
	c.2 = corpus callosum
	d.3 = hippocampus
	e.4 = thalamus 
	f.5 = putamen
	g.0 = out of brain
	h.White matter = regions 2, 4, 5
	i.Gray Matter = regions 1 & 3
