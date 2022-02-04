from pyteomics import mgf
import sys

### example of usage: python parse_mzmine_mgfs.py ~/Dropbox/workshop-FEB2022/mzmine_files/PoDP_div_genomes-TFL210409.mgf \
# ~/Dropbox/workshop-FEB2022/mzmine_files/ ~/Dropbox/workshop-FEB2022/NPOmix_python/gnps_names_list-220116.txt 300

input_file = sys.argv[1]
output_folder = sys.argv[2]
gnps_names_file= sys.argv[3]
mass_limit = sys.argv[4]

ccms_list = []
infile = open(gnps_names_file,'r')
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