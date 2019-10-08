MasterFunction<-function(tree_file, target_taxa, Fasta_path){
  
  fasta_files<-list.files(Fasta_path)
  
  
  for (i in 1:length(fasta_files)){
    ImportedTree<-read_annotated(tree_file, format="newick")
    
    MasterTaxa<-read.csv(target_taxa, header=FALSE)
    MasterTaxaL<-MasterTaxa$V1
    MasterTaxaList<-as.character(MasterTaxaL)
    GeneTaxa<-read.FASTA(file=paste(Fasta_path,fasta_files[1], sep=""))
    
    
    BigTreeTaxa<-ImportedTree$tip.label
    TaxaToDrop_2<-c()
    for (j in 1:length(BigTreeTaxa)) {
      if (BigTreeTaxa[j] %in% names(GeneTaxa)) {}
      else {TaxaToDrop_2<-list.append(TaxaToDrop_2, BigTreeTaxa[j])}
    }
    PrunedTree<-drop.tip(ImportedTree, TaxaToDrop_2)
    PrunedTreeTaxa<-c(PrunedTree$tip.label)
    TaxaToAnnotate<-c()
    for (k in 1:length(MasterTaxaList)){
      if (MasterTaxaList[k] %in% PrunedTreeTaxa)
        TaxaToAnnotate<-list.append(TaxaToAnnotate,MasterTaxaList[k])
    }
    MRCAout<-getMRCA(PrunedTree, TaxaToAnnotate)
    x<-print(MRCAout)
    NodesLoopList<-getDescendants(PrunedTree, x, curr=NULL)
    for (l in (NodesLoopList)){
      PrunedTree$node.comment[l]<-"{test}"
    }
    write_annotated(PrunedTree, file=paste(fasta_files[i],".nex", sep=""), format="nexus")
    print(PrunedTree)
  }
}
MasterFunction(tree_file="~/Academics/TAMUCC/HyPhy/TestMasterFunctionData/RAxML_bestTree.GobioEdit1_output", target_taxa = "~/Academics/TAMUCC/HyPhy/TestMasterFunctionData/MasterTaxaList.csv", Fasta_path = "~/Academics/TAMUCC/HyPhy/TestMasterFunctionData/Alignments3/")  
PrunedTree
#For some reason there are now 2 files appearing in the "Output" folder, with the names being "alignmentfileitwasmadefrom.fna.nwk
##NEVERMIND THAT'S BECAUSE THERE'S 2 FILES IN THE INPUT ALIGNMENT FOLDER AND THIS DOES IT TO EVERYTHING IN THE INPUT ALIGNMENT FOLDER
###The outputted files from this are named "alignmentfileitwasmadefrom.nwk" ...so same naming convention as it ever was.
#### But now they are nexus files with a taxa list, some sort of "trees" section that has numbers before the taxa names, then a newick tree at the very end
#####HYPHY/RELAX WORKS EVEN LESS ON THIS
