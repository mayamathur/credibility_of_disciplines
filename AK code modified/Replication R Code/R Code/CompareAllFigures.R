# Compare plots


setwd("/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication/FiguresandTables")
Filenames<-list.files(pattern="*.pdf")
Filenames2<-list.files(path="/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication/FiguresandTables",
  pattern="*.pdf")
sink(paste0("/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication/FiguresandTables/",  "CompareAllFigures.tex"));
for (name in intersect(Filenames,Filenames2)) {
    cat(sprintf("\\begin{figure} \n"));
    cat(sprintf("\\begin{subfigure}[b]{0.4\\textwidth} \n"));
    name1<-paste0(" \\includegraphics[width=\\textwidth]{",
                  "/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication/FiguresandTables/",name,
                  "} \n")
    name2<-paste0(" \\includegraphics[width=\\textwidth]{",
                  "/Users/virasemenora/Dropbox (MIT)/Isaiah Andrews Replication/FiguresandTables/",
                  name,
                  "} \n")
    cat(sprintf(name1));
    cat(sprintf("\\caption{Vira} \n"));
    cat(sprintf("\\label{fig:Vira} \n"));
    cat(sprintf("\\end{subfigure} \n"));
    cat(sprintf("\\begin{subfigure}[b]{0.4\\textwidth} \n"));
    cat(sprintf(name2));
    cat(sprintf("\\caption{Original} \n"));
    cat(sprintf("\\label{fig:Original} \n"));
    cat(sprintf("\\end{subfigure} \n"));
    cat(sprintf("\\end{figure} \n"));

}
#sink();
