# Author: Zhengwu Zhang

library(R.matlab)
library(glmnet)

rm(list=ls())
getwd();
#load data
# Zhengwu's Dropbox
setwd('Tensor_Decomposition/data')

# read csv file
table1<-read.csv("RESTRICTED_sony5650_6_18_2017_10_24_27.csv",header = TRUE);
table2<-read.csv("unrestricted_sony5650_6_18_2017_10_23_3.csv",header = TRUE);

#number of subject
# table1 - restricted
# table2 - unrestricted
nTol = length(table2$Subject)
AllSubject = table2$Subject
processed_sub_ids<-readMat('HCP1200ID.mat');
network_id = processed_sub_ids$hcp1200yesdiffusionID


####################################################################################
# choose the measures that related to the congnitive function of human brain
####################################################################################

allsubject_id = AllSubject;
handness = table1$Handedness;

#cog measure - response
cog_measure = matrix(0,175,length(AllSubject));

# confounding variables - 
confond_variable = matrix(0,2,length(AllSubject));

############################# Some Confonding Variables #######################
confond_variable[1,] = table1$Age_in_Yrs;
confond_variable[2,] = table2$Gender;


############################# Cognition related #######################
#Fluid Intelligence (Penn Progressive Matrices)
cog_measure[1,] = table2$PMAT24_A_CR 
cog_measure[2,] = table2$PMAT24_A_SI 
cog_measure[3,] = table2$PMAT24_A_RTCR #normal dist

#Instrument: Language/Reading Decoding (Oral Reading Recognition)
cog_measure[4,] = table2$ReadEng_Unadj #normal dist
cog_measure[5,] = table2$ReadEng_AgeAdj #normal dist

#Language/Vocabulary Comprehension (Picture Vocabulary)
cog_measure[6,] = table2$PicVocab_Unadj #normal dist
cog_measure[7,] = table2$PicVocab_AgeAdj #normal dist

#Instrument: Verbal Episodic Memory (Penn Word Memory Test)
cog_measure[8,] = table2$IWRD_TOT
cog_measure[9,] = table2$IWRD_RTC

#Processing Speed (Pattern Completion Processing Speed)
cog_measure[10,] = table2$ProcSpeed_Unadj #normal dist
cog_measure[11,] = table2$ProcSpeed_AgeAdj #normal dist
 
# Delay Discounting
cog_measure[12,] = table2$DDisc_SV_1mo_200 
cog_measure[13,] = table2$DDisc_SV_6mo_200 
cog_measure[14,] = table2$DDisc_SV_1yr_200 
cog_measure[15,] = table2$DDisc_SV_3yr_200 
cog_measure[16,] = table2$DDisc_SV_5yr_200 
cog_measure[17,] = table2$DDisc_SV_10yr_200 
cog_measure[18,] = table2$DDisc_SV_6mo_40K 
cog_measure[19,] = table2$DDisc_SV_1yr_40K 
cog_measure[20,] = table2$DDisc_SV_3yr_40K 
cog_measure[21,] = table2$DDisc_SV_5yr_40K 
cog_measure[22,] = table2$DDisc_SV_10yr_40K 
cog_measure[23,] = table2$DDisc_AUC_200 
cog_measure[24,] = table2$DDisc_AUC_40K #normal dist

#Instrument: Spatial Orientation (Variable Short Penn Line Orientation Test)
cog_measure[25,] = table2$VSPLOT_TC  #normal dist
cog_measure[26,] = table2$VSPLOT_CRTE #normal dist
cog_measure[27,] = table2$VSPLOT_OFF 

#Instrument: Sustained Attention (Short Penn Continuous Performance Test)
cog_measure[28,] = table2$SCPT_TP 
cog_measure[29,] = table2$SCPT_TN 
cog_measure[30,] = table2$SCPT_FP 
cog_measure[31,] = table2$SCPT_FN 
cog_measure[32,] = table2$SCPT_TPRT  #normal dist
cog_measure[33,] = table2$SCPT_SEN  
cog_measure[34,] = table2$SCPT_SPEC 
cog_measure[35,] = table2$SCPT_LRNR 

# Working Memory (List Sorting)
cog_measure[36,] = table2$ListSort_AgeAdj #normal dist
cog_measure[37,] = table2$ListSort_Unadj

#  Episodic Memory (Picture Sequence Memory) - C:Cognition
cog_measure[38,] = table2$PicSeq_AgeAdj #normal dist
cog_measure[39,] = table2$PicSeq_Unadj #normal dist

#Education and Income
cog_measure[40,] = table1$SSAGA_Educ
cog_measure[41,] = table1$SSAGA_Income

#Instrument: Executive Function/Cognitive Flexibility (Dimensional Change Card Sort)
cog_measure[42,] = table2$CardSort_Unadj #normal dist
cog_measure[43,] = table2$CardSort_AgeAdj

#Instrument: Executive Function/Inhibition (Flanker Task)
cog_measure[44,] = table2$Flanker_Unadj #normal dist
cog_measure[45,] = table2$Flanker_AgeAdj #normal dist


############################# Motor related #######################
# Instrument: Dexterity (9-hole Pegboard) - C: Motor
cog_measure[46,] = table2$Dexterity_AgeAdj  #normal dist
cog_measure[47,] = table2$Dexterity_Unadj  #normal dist

# Instrument: Endurance (2 minute walk test)
cog_measure[48,] = table2$Endurance_Unadj  #normal dist
cog_measure[49,] = table2$Endurance_AgeAdj  #normal dist

# Instrument: Locomotion (4-meter walk test)
cog_measure[50,] = table2$GaitSpeed_Comp  #normal dist

#Instrument: Strength (Grip Strength Dynamometry)
cog_measure[51,] = table2$Strength_Unadj  #normal dist
cog_measure[52,] = table2$Strength_AgeAdj #normal dist



############################# Substance Use #######################
#Instrument: Alcohol Use and Dependence
cog_measure[53,] = table1$SSAGA_Alc_D4_Dp_Sx 
cog_measure[54,] = table1$SSAGA_Alc_D4_Ab_Dx 
cog_measure[55,] = table1$SSAGA_Alc_D4_Ab_Sx 
cog_measure[56,] = table1$SSAGA_Alc_D4_Dp_Dx 
cog_measure[57,] = table1$SSAGA_Alc_12_Drinks_Per_Day 
cog_measure[58,] = table1$SSAGA_Alc_12_Frq
cog_measure[59,] = table1$SSAGA_Alc_12_Frq_5plus
cog_measure[60,] = table1$SSAGA_Alc_12_Frq_Drk
cog_measure[61,] = table1$SSAGA_Alc_12_Max_Drinks
cog_measure[62,] = table1$SSAGA_Alc_Age_1st_Use
cog_measure[63,] = table1$SSAGA_Alc_Hvy_Drinks_Per_Day
cog_measure[64,] = table1$SSAGA_Alc_Hvy_Frq
cog_measure[65,] = table1$SSAGA_Alc_Hvy_Frq_5plus
cog_measure[66,] = table1$SSAGA_Alc_Hvy_Frq_Drk
cog_measure[67,] = table1$SSAGA_Alc_Hvy_Max_Drinks

#Instrument: Tobacco Use and Dependence
cog_measure[68,] = table1$SSAGA_TB_Age_1st_Cig 
cog_measure[69,] = table1$SSAGA_TB_DSM_Difficulty_Quitting 
cog_measure[70,] = table1$SSAGA_TB_Max_Cigs 
cog_measure[71,] = table1$SSAGA_TB_Reg_CPD 
cog_measure[72,] = table1$SSAGA_TB_Yrs_Smoked 
cog_measure[73,] = table1$SSAGA_TB_Still_Smoking
cog_measure[74,] = table1$Times_Used_Any_Tobacco_Today
cog_measure[75,] = table1$Avg_Weekend_Any_Tobacco_7days
cog_measure[76,] = table1$Total_Cigars_7days
cog_measure[77,] = table1$Avg_Weekend_Cigars_7days
cog_measure[78,] = table1$Num_Days_Used_Any_Tobacco_7days

#Instrument: Illicit Drug Use
cog_measure[79,] = table1$SSAGA_Times_Used_Illicits 
cog_measure[80,] = table1$SSAGA_Times_Used_Cocaine 
cog_measure[81,] = table1$SSAGA_Times_Used_Hallucinogens 
cog_measure[82,] = table1$SSAGA_Times_Used_Opiates 
cog_measure[83,] = table1$SSAGA_Times_Used_Sedatives 
cog_measure[84,] = table1$SSAGA_Times_Used_Stimulants

#Instrument: Marijuana Use and Dependence
cog_measure[85,] = table1$SSAGA_Mj_Use 
cog_measure[86,] = table1$SSAGA_Mj_Ab_Dep 
cog_measure[87,] = table1$SSAGA_Mj_Age_1st_Use 
cog_measure[88,] = table1$SSAGA_Mj_Times_Used 

############################# Psychiatric and Life Function #######################
# lift function
cog_measure[89,] = table1$ASR_Anxd_Raw 
cog_measure[90,] = table1$ASR_Anxd_Pct 
cog_measure[91,] = table1$ASR_Witd_Raw 
cog_measure[92,] = table1$ASR_Witd_Pct 
cog_measure[93,] = table1$ASR_Soma_Raw 
cog_measure[94,] = table1$ASR_Soma_Pct
cog_measure[95,] = table1$ASR_Thot_Raw 
cog_measure[96,] = table1$ASR_Thot_Pct
cog_measure[97,] = table1$ASR_Attn_Raw 
cog_measure[98,] = table1$ASR_Attn_Pct
cog_measure[99,] = table1$ASR_Aggr_Raw 
cog_measure[100,] = table1$ASR_Aggr_Pct
cog_measure[101,] = table1$ASR_Rule_Raw 
cog_measure[102,] = table1$ASR_Rule_Pct
cog_measure[103,] = table1$ASR_Intr_Raw 
cog_measure[104,] = table1$ASR_Intr_Pct #normal dist
cog_measure[105,] = table1$ASR_Oth_Raw 
cog_measure[106,] = table1$ASR_Crit_Raw
cog_measure[107,] = table1$ASR_Intn_Raw 
cog_measure[108,] = table1$ASR_Intn_T #normal dist
cog_measure[109,] = table1$ASR_Extn_Raw 
cog_measure[110,] = table1$ASR_Extn_T #normal dist
cog_measure[111,] = table1$ASR_TAO_Sum  #normal dist
cog_measure[112,] = table1$ASR_Totp_Raw #normal dist
cog_measure[113,] = table1$ASR_Totp_T #normal dist
cog_measure[114,] = table1$DSM_Depr_Raw 
cog_measure[115,] = table1$DSM_Depr_Pct
cog_measure[116,] = table1$DSM_Anxi_Raw
cog_measure[117,] = table1$DSM_Anxi_Pct
cog_measure[118,] = table1$DSM_Somp_Raw
cog_measure[119,] = table1$DSM_Somp_Pct
cog_measure[120,] = table1$DSM_Avoid_Raw
cog_measure[121,] = table1$DSM_Avoid_Pct 
cog_measure[122,] = table1$DSM_Adh_Raw
cog_measure[123,] = table1$DSM_Adh_Pct
cog_measure[124,] = table1$DSM_Hype_Raw
cog_measure[125,] = table1$DSM_Antis_Raw
cog_measure[126,] = table1$DSM_Antis_Pct

#Restricted Instrument: Psychiatric History
cog_measure[127,] = table1$SSAGA_ChildhoodConduct
cog_measure[128,] = table1$SSAGA_PanicDisorder
cog_measure[129,] = table1$SSAGA_Agoraphobia
cog_measure[130,] = table1$SSAGA_Depressive_Ep
cog_measure[131,] = table1$SSAGA_Depressive_Sx

############################# Sensory #######################
# Instrument: Audition (Words in Noise
cog_measure[132,] = table2$Noise_Comp 
cog_measure[133,] = table2$Odor_Unadj #Odor
cog_measure[134,] = table2$Odor_AgeAdj

#cog_measure[65,] = table1$PainIntens_RawScore #Odor
cog_measure[135,] = table2$PainInterf_Tscore

cog_measure[136,] = table2$Taste_Unadj #Taste #normal dist
cog_measure[137,] = table2$Taste_AgeAdj #normal dist

cog_measure[138,] = table1$EVA_Num #Vision
cog_measure[139,] = table1$EVA_Denom

#Instrument: Contrast Sensitivity (Mars Contrast Sensitivity)
cog_measure[140,] = table2$Mars_Log_Score #Vision
cog_measure[141,] = table2$Mars_Errs
cog_measure[142,] = table2$Mars_Final

############################# Emotion #######################
# Emotion recognition
cog_measure[143,] = table2$ER40_CR 
cog_measure[144,] = table2$ER40_CRT #normal dist
cog_measure[145,] = table2$ER40ANG 
cog_measure[146,] = table2$ER40FEAR
cog_measure[147,] = table2$ER40HAP 
cog_measure[148,] = table2$ER40NOE
cog_measure[149,] = table2$ER40SAD

# Instrument: Negative Affect 
cog_measure[150,] = table2$AngAffect_Unadj #normal dist
cog_measure[151,] = table2$AngHostil_Unadj  #normal dist
cog_measure[152,] = table2$AngAggr_Unadj
cog_measure[153,] = table2$FearAffect_Unadj #normal dist
cog_measure[154,] = table2$FearSomat_Unadj
cog_measure[155,] = table2$Sadness_Unadj #normal dist

# Instrument: Stress and Self Efficacy
cog_measure[156,] = table2$PercStress_Unadj #normal dist
cog_measure[157,] = table2$SelfEff_Unadj #normal dist

# Psychological Well-being 
cog_measure[158,] = table2$LifeSatisf_Unadj #normal dist
cog_measure[159,] = table2$MeanPurp_Unadj #normal dist
cog_measure[160,] = table2$PosAffect_Unadj

# Social Relationships 
cog_measure[161,] = table2$Friendship_Unadj #normal dist
cog_measure[162,] = table2$Loneliness_Unadj
cog_measure[163,] = table2$PercHostil_Unadj #normal dist
cog_measure[164,] = table2$PercReject_Unadj
cog_measure[165,] = table2$EmotSupp_Unadj
cog_measure[166,] = table2$InstruSupp_Unadj

############################# Personality #######################
# Instrument: Five Factor Model (NEO-FFI)
cog_measure[167,] = table2$NEOFAC_A  #normal dist
cog_measure[168,] = table2$NEOFAC_O  #normal dist
cog_measure[169,] = table2$NEOFAC_C  #normal dist
cog_measure[170,] = table2$NEOFAC_N #normal dist
cog_measure[171,] = table2$NEOFAC_E  #normal dist


############################# Health and Family History #######################
cog_measure[172,] = table1$Height #normal dist
cog_measure[173,] = table1$Weight  #normal dist
cog_measure[174,] = table1$BMI #normal dist
cog_measure[175,] = table2$PSQI_Score



writeMat('data/HCP_Covariates.mat',confond_variable = confond_variable, cog_measure = cog_measure,allsubject_id = allsubject_id,handness = handness)
