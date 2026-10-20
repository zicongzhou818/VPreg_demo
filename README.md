# VPreg

VPreg is a non-rigid image registration method, it is required to have your samples preprocessed. The protocols are included and ordered as follows:
1. N4 bias correction;
2. Skull skripping or brain extraction;
3. Delineation, i.e., linear registration from samples to template. 

To see this testing example, download the codes and files and update your working path in test_SymVPcontrol_sample2template.m, before running it.
This package gives a test on registrations from samples to a template. It produces the forward and and backward mappings, as well as their corresponding Mapped and inverseMapped images. 
