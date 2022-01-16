from pyteomics import mgf
import sys

input_file = sys.argv[1]
output_folder = sys.argv[2]
mass_limit = sys.argv[3]

ccms_list = []
infile = open('/Users/tiagoferreiraleao/Dropbox/tiago-NAS/NPOmix_python/gnps_names_list-220116.txt','r')
for line in infile:
    ccms_list.append(line.rstrip('\n'))
    
count = 0

with mgf.MGF(input_file) as reader:
    for spectrum in reader:
        for key in spectrum.keys():
            if key == 'params':
                if spectrum[key]['pepmass'][0] > float(mass_limit):
                    count += 1
                    individual_spec = []
                    individual_spec.append(spectrum)
                    output_file = '%s%s.mgf'%(output_folder,ccms_list[count])
                    mgf.write(spectra=individual_spec, header='', output=output_file)

print ('%s files extracted and renamed'%count)