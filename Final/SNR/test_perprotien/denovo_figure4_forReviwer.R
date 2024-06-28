
data=read.csv("C:\\Workplace\\Python\\AnalysisForThePaper\\NEH\\Final\\Simulation_ape_mpe_asymp\\_data_mpe_ape_3.csv")

pxt_exp=c(rep(-1, dim(data)[1] ))
neh_exp=c(rep(-1, dim(data)[1] ))
pxt_theo=c(rep(-1, dim(data)[1] ))
neh_theo=c(rep(-1, dim(data)[1] ))

for (i in seq(1,(dim(data)[1]))){
  
  vals_exp=Compute_neh_pxt(data$DataRecord.I0[i],
                           data$DataRecord.I1[i],
                           data$DataRecord.I2[i],                           
                           data$DataRecord.I0_exp_t[i],
                           data$DataRecord.I1_exp_t[i],
                           data$DataRecord.I2_exp_t[i])
  
  vals_theo= c(0,0) 
  
  # print(c("px(t)",vals[1],"NEH",vals[2]))
  pxt_theo[i]=vals_theo[1];
  neh_theo[i]=vals_theo[2]
  
  pxt_exp[i]=vals_exp[1];
  neh_exp[i]=vals_exp[2]
  
}

df=data.frame(data$DataRecord.NEH,c(rep(0.046,dim(data)[1])),pxt_theo,neh_theo,pxt_exp,neh_exp)
colnames(df)=c('T_NEH','T_PXT','PXT_THEO','NEH_THEO','PXT_EXP','NEH_EXP')
write.csv(df,"denovo_figure4.csv",row.names = FALSE)
