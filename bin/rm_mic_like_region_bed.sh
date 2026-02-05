
bed="GCA_027405865.1_94003_Monoecious_Salix_purpurea_Haplotype_2_assembly.bed"

bam="/cfs/klemming/projects/snic/snic2022-6-144/CHENYU/chap3_popgen/251028/bam_nuclear/ERR10493277.rm_low_complex.GCA_027405865.1_94003_Monoecious_Salix_purpurea_Haplotype_2_assembly.bam"

out="${bam%.bam}.mic_masked.bam"
    
#output the bam file that is not in the regions
samtools view -b -h -L $bed_dir/"$bed" -U $out -o /dev/null "$bam" 
