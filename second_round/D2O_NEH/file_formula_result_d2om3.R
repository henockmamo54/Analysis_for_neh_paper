
data=read.csv("C:\\Workplace\\C++\\d2ome_restructure\\d2ome_restructure\\d2ome_GUI\\d2ome_GUI\\v2\\bin\\Debug\\Estimated_neh_numbers_utmb.csv")

pxt=c(rep(-1, dim(data)[1] ))
neh=c(rep(-1, dim(data)[1] ))

for (i in seq(1,(dim(data)[1]))){
  
  vals=Compute_neh_pxt(data$i0_0[i],data$i1_0[i],data$i2_0[i],
                       data$i0_31[i],data$i1_31[i],data$i2_31[i]); # computes the px(t) value and 
  #the corresponding NEH value
  
  print(c("px(t)",vals[1],"NEH",vals[2]))
  pxt[i]=vals[1];
  neh[i]=vals[2]
  
}

df <- data.frame(data$Peptide,data$Protein,data$Charge,data$T_NEH,pxt,neh)
write.csv(df, "neh_from_formula_utmb.csv", row.names = FALSE)
print('done!')
