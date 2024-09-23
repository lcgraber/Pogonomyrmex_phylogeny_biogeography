##sortadate code used for selection of clock-like loci##

#https://github.com/FePhyFoFum/SortaDate?tab=readme-ov-file
# The expected input is a set of genes that have been aligned and for which there are gene trees

export PATH="/home/lcg65/Desktop/phyx-master/phyx-master/src:$PATH"

#I think I have to root trees before doing any of these

#rooting gene trees
for i in /home/lcg65/Desktop/PogoUCE/aug_2_time_tree_internal_trimmed_2024_gblocks_clean_75p/all_tree_for_astral/*.treefile; 
do 
    name=$(basename "$i" .treefile)
    pxrr -t $i -g hylomyrma_blandiens_VEN_AG_EX1598 -o /home/lcg65/Desktop/PogoUCE/SortaDate_rootedtrees/$name ;
done

#rooting species tree
pxrr -t /home/lcg65/Desktop/PogoUCE/ASTRAL_runs/20August2024/all_gene_time_tree.treefile -g hylomyrma_blandiens_VEN_AG_EX1598 -o /home/lcg65/Desktop/PogoUCE/SortaDate_rootedtrees/speciestree.treefile

#get root to tip variance
python /home/lcg65/Desktop/SortaDate-master/src/get_var_length.py /home/lcg65/Desktop/PogoUCE/SortaDate_rootedtrees --flend .nexus --outf variance --outg hylomyrma_blandiens_VEN_AG_EX1598

#get bipartition support
python /home/lcg65/Desktop/SortaDate-master/src/get_bp_genetrees.py /home/lcg65/Desktop/PogoUCE/SortaDate_rootedtrees/ /home/lcg65/Desktop/PogoUCE/SortaDate_rootedtrees/species_tree/speciestree.treefile --flend .nexus --outf bp_support 

#combine results
python /home/lcg65/Desktop/SortaDate-master/src/combine_results.py variance bp_support --outf comb

#sort adn get list of good genes
python /home/lcg65/Desktop/SortaDate-master/src/get_good_genes.py comb --max 100 --order 3,1,2 --outf gg
