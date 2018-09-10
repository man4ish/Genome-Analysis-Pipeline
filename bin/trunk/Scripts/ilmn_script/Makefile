CXX = g++
CFLAGS=-Wall -O3 -Werror
LFLAGS=-lpthread

OBJ_DIR = obj
OUT_DIR = target

HOME_DIR = ../..

TARGET_SXReadSNV = $(OUT_DIR)/SXReadSNV_Ilm
SRC_DIR_SXReadSNV = $(HOME_DIR)/SXReadSNV_Ilm

TARGET_SSNpParser = $(OUT_DIR)/SSNpParser
SRC_DIR_SSNpParser  = $(HOME_DIR)/SSNpParser

TARGET_SXFilterSRB4Merge =  $(OUT_DIR)/SXFilterSRB4Merge_SingleSample
SRC_DIR_SXFilterSRB4Merge = $(HOME_DIR)/SXFilterSRB4Merge

TARGET_SXPrintZygosityCNV = $(OUT_DIR)/SXPZygoCNV_SingleSampleEx
TARGET_SXPZygoINDEL = $(OUT_DIR)/SXPZygoINDEL_SingleSample
TARGET_SXPrintCNV = $(OUT_DIR)/SXPrintCNVNS
SRC_DIR_SXPrintZygosityCNV = $(HOME_DIR)/SXPrintZygosityCNV

TARGET_SXAppendScoreSNP = $(OUT_DIR)/SXAppendScoreSNP_Ilmn
TARGET_SXAppendScoreINS = $(OUT_DIR)/SXAppendScoreINS_Ilmn
TARGET_SXAppendScoreDEL = $(OUT_DIR)/SXAppendScoreDEL_Ilmn
TARGET_SXAppendScoreSNP_Ilmn_B4Annotation = $(OUT_DIR)/SXAppendScoreSNP_Ilmn_B4Annotation
TARGET_SXAppendScoreINS_Ilmn_B4Annotation = $(OUT_DIR)/SXAppendScoreINS_Ilmn_B4Annotation
TARGET_SXAppendScoreDEL_Ilmn_B4Annotation = $(OUT_DIR)/SXAppendScoreDEL_Ilmn_B4Annotation
SRC_DIR_SXAppendScore = $(HOME_DIR)/SXAppendScore

TARGET_SXAppendDrugInfo = $(OUT_DIR)/SXAppendDrugInfo
SRC_DIR_SXAppendDrugInfo = $(HOME_DIR)/SXAppendDrugInfo

TARGET_SXAppendPerfectDenPoly = $(OUT_DIR)/SXAppendPDenMinus2Poly_SingleSample
SRC_DIR_SXAppendPerfectDenPoly = $(HOME_DIR)/SXAppendPerfectDenPoly
                                             
TARGET_SXSynoMapEx = $(OUT_DIR)/SXSynoMapEx
SRC_DIR_SXSynoMapEx = $(HOME_DIR)/SXSynoMapEx

TARGET_Avg_Confi_Score = $(OUT_DIR)/Avg_Confi_Score
SRC_DIR_Avg_Confi_Score = $(HOME_DIR)/Avg_Confi_Score

TARGET_SXFilterByCScore =  $(OUT_DIR)/SXFilterByCScore
SRC_DIR_SXFilterByCScore = $(HOME_DIR)/SXFilterByCScore

TARGET_SXPrintStatsSNP = $(OUT_DIR)/SXPrintStatsSNP
TARGET_SXPrintStatsINS = $(OUT_DIR)/SXPrintStatsINS
TARGET_SXPrintStatsDEL = $(OUT_DIR)/SXPrintStatsDEL
SRC_DIR_SXPrintStats = $(HOME_DIR)/SXPrintStats

all: targetdirs $(TARGET_SXReadSNV) $(TARGET_SXFilterSRB4Merge) $(TARGET_SXPrintZygosityCNV) \
	$(TARGET_SXAppendScoreSNP) $(TARGET_SXAppendScoreINS) $(TARGET_SXAppendScoreDEL) \
	$(TARGET_SXAppendScoreSNP_Ilmn_B4Annotation) $(TARGET_SXAppendScoreINS_Ilmn_B4Annotation) \
	$(TARGET_SXAppendScoreDEL_Ilmn_B4Annotation) $(TARGET_SXAppendDrugInfo) \
	$(TARGET_SSNpParser) $(TARGET_SXAppendPerfectDenPoly) $(TARGET_SXPZygoINDEL) \
	$(TARGET_SXPrintCNV) $(TARGET_SXSynoMapEx) $(TARGET_Avg_Confi_Score) \
	$(TARGET_SXFilterByCScore) $(TARGET_SXPrintStatsSNP) $(TARGET_SXPrintStatsINS) $(TARGET_SXPrintStatsDEL)

#----------------------------------------------------------------------------------------------------------------------------
$(TARGET_SXReadSNV): $(OBJ_DIR)/main.o $(OBJ_DIR)/CSXReadSNV.o $(OBJ_DIR)/SXDBSNPgChecker.o $(OBJ_DIR)/SXDBIndelgChecker.o $(OBJ_DIR)/SXAnnotate.o $(OBJ_DIR)/SXDBIndelgChecker_OdbSnpfmt.o
	$(CXX) $(CFLAGS) $(LFLAGS) -o $(TARGET_SXReadSNV) $(OBJ_DIR)/main.o $(OBJ_DIR)/CSXReadSNV.o $(OBJ_DIR)/SXDBSNPgChecker.o $(OBJ_DIR)/SXDBIndelgChecker.o $(OBJ_DIR)/SXAnnotate.o $(OBJ_DIR)/SXDBIndelgChecker_OdbSnpfmt.o

$(OBJ_DIR)/main.o: $(SRC_DIR_SXReadSNV)/main.cpp $(SRC_DIR_SXReadSNV)/CSXReadSNV.h
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/main.o $(SRC_DIR_SXReadSNV)/main.cpp

$(OBJ_DIR)/CSXReadSNV.o: $(SRC_DIR_SXReadSNV)/CSXReadSNV.cpp $(SRC_DIR_SXReadSNV)/CSXReadSNV.h $(SRC_DIR_SXReadSNV)/typedef.h $(SRC_DIR_SXReadSNV)/SXDBSNPgChecker.h $(SRC_DIR_SXReadSNV)/SXDBIndelgChecker.h $(SRC_DIR_SXReadSNV)/SXAnnotate.h	
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/CSXReadSNV.o $(SRC_DIR_SXReadSNV)/CSXReadSNV.cpp

$(OBJ_DIR)/SXDBSNPgChecker.o: $(SRC_DIR_SXReadSNV)/SXDBSNPgChecker.cpp $(SRC_DIR_SXReadSNV)/SXDBSNPgChecker.h
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXDBSNPgChecker.o $(SRC_DIR_SXReadSNV)/SXDBSNPgChecker.cpp

$(OBJ_DIR)/SXDBIndelgChecker.o: $(SRC_DIR_SXReadSNV)/SXDBIndelgChecker.cpp $(SRC_DIR_SXReadSNV)/SXDBIndelgChecker.h
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXDBIndelgChecker.o $(SRC_DIR_SXReadSNV)/SXDBIndelgChecker.cpp

$(OBJ_DIR)/SXAnnotate.o: $(SRC_DIR_SXReadSNV)/SXAnnotate.cpp $(SRC_DIR_SXReadSNV)/SXAnnotate.h
	${CXX} $(CFLAGS) -c -o $(OBJ_DIR)/SXAnnotate.o  $(SRC_DIR_SXReadSNV)/SXAnnotate.cpp

$(OBJ_DIR)/SXDBIndelgChecker_OdbSnpfmt.o: $(SRC_DIR_SXReadSNV)/SXDBIndelgChecker_OdbSnpfmt.cpp $(SRC_DIR_SXReadSNV)/SXDBIndelgChecker_OdbSnpfmt.h
	${CXX} $(CFLAGS) -c -o $(OBJ_DIR)/SXDBIndelgChecker_OdbSnpfmt.o $(SRC_DIR_SXReadSNV)/SXDBIndelgChecker_OdbSnpfmt.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SSNpParser): $(OBJ_DIR)/SSNpParser.o
	$(CXX) $(CFLAGS) -o $(TARGET_SSNpParser) $(OBJ_DIR)/SSNpParser.o

$(OBJ_DIR)/SSNpParser.o: $(SRC_DIR_SSNpParser)/SSNpParser.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SSNpParser.o $(SRC_DIR_SSNpParser)/SSNpParser.cpp

#----------------------------------------------------------------------------------------------------------------------------   

$(TARGET_SXAppendPerfectDenPoly): $(OBJ_DIR)/SXAppendPDenMinus2Poly_SingleSample.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXAppendPerfectDenPoly) $(OBJ_DIR)/SXAppendPDenMinus2Poly_SingleSample.o  

$(OBJ_DIR)/SXAppendPDenMinus2Poly_SingleSample.o: $(SRC_DIR_SXAppendPerfectDenPoly)/SXAppendPDenMinus2Poly_SingleSample.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXAppendPDenMinus2Poly_SingleSample.o $(SRC_DIR_SXAppendPerfectDenPoly)/SXAppendPDenMinus2Poly_SingleSample.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXFilterSRB4Merge): $(OBJ_DIR)/SXFilterSRB4Merge_SingleSample.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXFilterSRB4Merge) $(OBJ_DIR)/SXFilterSRB4Merge_SingleSample.o

$(OBJ_DIR)/SXFilterSRB4Merge_SingleSample.o: $(SRC_DIR_SXFilterSRB4Merge)/SXFilterSRB4Merge_SingleSample.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXFilterSRB4Merge_SingleSample.o $(SRC_DIR_SXFilterSRB4Merge)/SXFilterSRB4Merge_SingleSample.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXPrintZygosityCNV): $(OBJ_DIR)/SXPZygoCNV_SingleSampleEx.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXPrintZygosityCNV) $(OBJ_DIR)/SXPZygoCNV_SingleSampleEx.o

$(OBJ_DIR)/SXPZygoCNV_SingleSampleEx.o: $(SRC_DIR_SXPrintZygosityCNV)/SXPZygoCNV_SingleSampleEx.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXPZygoCNV_SingleSampleEx.o $(SRC_DIR_SXPrintZygosityCNV)/SXPZygoCNV_SingleSampleEx.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXAppendScoreSNP): $(OBJ_DIR)/SXAppendScoreSNP_Ilmn.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXAppendScoreSNP) $(OBJ_DIR)/SXAppendScoreSNP_Ilmn.o

$(OBJ_DIR)/SXAppendScoreSNP_Ilmn.o: $(SRC_DIR_SXAppendScore)/SXAppendScoreSNP_Ilmn.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXAppendScoreSNP_Ilmn.o $(SRC_DIR_SXAppendScore)/SXAppendScoreSNP_Ilmn.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXAppendScoreINS): $(OBJ_DIR)/SXAppendScoreINS_Ilmn.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXAppendScoreINS) $(OBJ_DIR)/SXAppendScoreINS_Ilmn.o

$(OBJ_DIR)/SXAppendScoreINS_Ilmn.o: $(SRC_DIR_SXAppendScore)/SXAppendScoreINS_Ilmn.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXAppendScoreINS_Ilmn.o $(SRC_DIR_SXAppendScore)/SXAppendScoreINS_Ilmn.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXAppendScoreDEL): $(OBJ_DIR)/SXAppendScoreDEL_Ilmn.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXAppendScoreDEL) $(OBJ_DIR)/SXAppendScoreDEL_Ilmn.o

$(OBJ_DIR)/SXAppendScoreDEL_Ilmn.o: $(SRC_DIR_SXAppendScore)/SXAppendScoreDEL_Ilmn.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXAppendScoreDEL_Ilmn.o $(SRC_DIR_SXAppendScore)/SXAppendScoreDEL_Ilmn.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXAppendScoreSNP_Ilmn_B4Annotation): $(OBJ_DIR)/SXAppendScoreSNP_Ilmn_B4Annotation.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXAppendScoreSNP_Ilmn_B4Annotation) $(OBJ_DIR)/SXAppendScoreSNP_Ilmn_B4Annotation.o

$(OBJ_DIR)/SXAppendScoreSNP_Ilmn_B4Annotation.o: $(SRC_DIR_SXAppendScore)/SXAppendScoreSNP_Ilmn_B4Annotation.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXAppendScoreSNP_Ilmn_B4Annotation.o $(SRC_DIR_SXAppendScore)/SXAppendScoreSNP_Ilmn_B4Annotation.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXAppendScoreINS_Ilmn_B4Annotation): $(OBJ_DIR)/SXAppendScoreINS_Ilmn_B4Annotation.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXAppendScoreINS_Ilmn_B4Annotation) $(OBJ_DIR)/SXAppendScoreINS_Ilmn_B4Annotation.o

$(OBJ_DIR)/SXAppendScoreINS_Ilmn_B4Annotation.o: $(SRC_DIR_SXAppendScore)/SXAppendScoreINS_Ilmn_B4Annotation.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXAppendScoreINS_Ilmn_B4Annotation.o $(SRC_DIR_SXAppendScore)/SXAppendScoreINS_Ilmn_B4Annotation.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXAppendScoreDEL_Ilmn_B4Annotation): $(OBJ_DIR)/SXAppendScoreDEL_Ilmn_B4Annotation.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXAppendScoreDEL_Ilmn_B4Annotation) $(OBJ_DIR)/SXAppendScoreDEL_Ilmn_B4Annotation.o

$(OBJ_DIR)/SXAppendScoreDEL_Ilmn_B4Annotation.o: $(SRC_DIR_SXAppendScore)/SXAppendScoreDEL_Ilmn_B4Annotation.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXAppendScoreDEL_Ilmn_B4Annotation.o $(SRC_DIR_SXAppendScore)/SXAppendScoreDEL_Ilmn_B4Annotation.cpp

#---------------------------------------------------------------------------------------------------------------------------
$(TARGET_SXAppendDrugInfo): $(OBJ_DIR)/SXAppendDrugInfo.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXAppendDrugInfo) $(OBJ_DIR)/SXAppendDrugInfo.o

$(OBJ_DIR)/SXAppendDrugInfo.o: $(SRC_DIR_SXAppendDrugInfo)/SXAppendDrugInfo.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXAppendDrugInfo.o $(SRC_DIR_SXAppendDrugInfo)/SXAppendDrugInfo.cpp

#----------------------------------------------------------------------------------------------------------------------------	

$(TARGET_SXPZygoINDEL): $(OBJ_DIR)/SXPZygoINDEL_SingleSample.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXPZygoINDEL) $(OBJ_DIR)/SXPZygoINDEL_SingleSample.o

$(OBJ_DIR)/SXPZygoINDEL_SingleSample.o: $(SRC_DIR_SXPrintZygosityCNV)/SXPZygoINDEL_SingleSample.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXPZygoINDEL_SingleSample.o $(SRC_DIR_SXPrintZygosityCNV)/SXPZygoINDEL_SingleSample.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXPrintCNV): $(OBJ_DIR)/SXPrintCNVNS.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXPrintCNV) $(OBJ_DIR)/SXPrintCNVNS.o

$(OBJ_DIR)/SXPrintCNVNS.o: $(SRC_DIR_SXPrintZygosityCNV)/SXPrintCNVNS.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXPrintCNVNS.o $(SRC_DIR_SXPrintZygosityCNV)/SXPrintCNVNS.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXSynoMapEx): $(OBJ_DIR)/SXSynoMapEx.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXSynoMapEx) $(OBJ_DIR)/SXSynoMapEx.o

$(OBJ_DIR)/SXSynoMapEx.o: $(SRC_DIR_SXSynoMapEx)/SXSynoMapEx.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXSynoMapEx.o $(SRC_DIR_SXSynoMapEx)/SXSynoMapEx.cpp

#----------------------------------------------------------------------------------------------------------------------------

$(TARGET_SXFilterByCScore): $(OBJ_DIR)/SXFilterByCScore.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXFilterByCScore) $(OBJ_DIR)/SXFilterByCScore.o

$(OBJ_DIR)/SXFilterByCScore.o: $(SRC_DIR_SXFilterByCScore)/SXFilterByCScore.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXFilterByCScore.o $(SRC_DIR_SXFilterByCScore)/SXFilterByCScore.cpp

#----------------------------------------------------------------------------------------------------------------------------
$(TARGET_Avg_Confi_Score): $(OBJ_DIR)/Avg_Confi_Score.o
	$(CXX) $(CFLAGS) -o $(TARGET_Avg_Confi_Score) $(OBJ_DIR)/Avg_Confi_Score.o

$(OBJ_DIR)/Avg_Confi_Score.o: $(SRC_DIR_Avg_Confi_Score)/Avg_Confi_Score.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/Avg_Confi_Score.o $(SRC_DIR_Avg_Confi_Score)/Avg_Confi_Score.cpp

#----------------------------------------------------------------------------------------------------------------------------
$(TARGET_SXPrintStatsSNP): $(OBJ_DIR)/SXPrintStatsSNP.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXPrintStatsSNP) $(OBJ_DIR)/SXPrintStatsSNP.o

$(OBJ_DIR)/SXPrintStatsSNP.o: $(SRC_DIR_SXPrintStats)/SXPrintStatsSNP.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXPrintStatsSNP.o $(SRC_DIR_SXPrintStats)/SXPrintStatsSNP.cpp

#----------------------------------------------------------------------------------------------------------------------------
$(TARGET_SXPrintStatsINS): $(OBJ_DIR)/SXPrintStatsINS.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXPrintStatsINS) $(OBJ_DIR)/SXPrintStatsINS.o

$(OBJ_DIR)/SXPrintStatsINS.o: $(SRC_DIR_SXPrintStats)/SXPrintStatsINS.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXPrintStatsINS.o $(SRC_DIR_SXPrintStats)/SXPrintStatsINS.cpp

#----------------------------------------------------------------------------------------------------------------------------
$(TARGET_SXPrintStatsDEL): $(OBJ_DIR)/SXPrintStatsDEL.o
	$(CXX) $(CFLAGS) -o $(TARGET_SXPrintStatsDEL) $(OBJ_DIR)/SXPrintStatsDEL.o

$(OBJ_DIR)/SXPrintStatsDEL.o: $(SRC_DIR_SXPrintStats)/SXPrintStatsDEL.cpp
	$(CXX) $(CFLAGS) -c -o $(OBJ_DIR)/SXPrintStatsDEL.o $(SRC_DIR_SXPrintStats)/SXPrintStatsDEL.cpp

#----------------------------------------------------------------------------------------------------------------------------
targetdirs:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(OUT_DIR)
