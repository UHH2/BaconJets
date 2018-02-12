#! /usr/bin/env python

import subprocess
import sys
import os
import os.path

output_name = "controlplots_defaultName.tex"
# print len(sys.argv)
# for a in sys.argv:
#     print a

if len(sys.argv) > 1:
    output_name = sys.argv[1] + ".tex"

img_path_base_2 = "/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17D_Data_woHOTVR_woXCON_test_JEC_V4/RunD_17Nov17_2017_newMC_woGen"
img_path_base_2 = img_path_base_2 + "/plots/"

if not os.path.isdir(img_path_base_2):
    raise RuntimeError("Path %s does not exist" % img_path_base_2)

img_path_base = img_path_base_2 + "control/"
img_path_base_D = img_path_base
img_path_baseMain_D = img_path_base_2

img_path_baseMain_E = "/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17E_Data_woHOTVR_woXCON_test_JEC_V1/RunE_17Nov17_2017_newMC_woGen/plots/"
img_path_base_E = img_path_baseMain_E + "control/"


img_path_base_th_D = "/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17D_Data_DeriveThresholds_woHOTVR_woXCONE/RunD_17Nov17_2017_Di/plots/thresholds"

img_path_base_th_E = "/nfs/dust/cms/user/garbersc/forBaconJets/17Nov2017/Residuals/Run17E_Data_DeriveThresholds_woHOTVR_woXCONE_JEC_V4/RunE_17Nov17_2017/plots/thresholds"

# img_path_base_B = "/nfs/dust/cms/user/garbersc/forBaconJets/2017PromptReco/Residuals/Run17BC_Data/RunB_5/plots/control/"
# img_path_base_C = "/nfs/dust/cms/user/garbersc/forBaconJets/2017PromptReco/Residuals/Run17BC_Data/RunC_5/plots/control/"
# img_path_baseMain_B = "/nfs/dust/cms/user/garbersc/forBaconJets/2017PromptReco/Residuals/Run17BC_Data/RunB_5/plots/"
# img_path_baseMain_C =
# "/nfs/dust/cms/user/garbersc/forBaconJets/2017PromptReco/Residuals/Run17BC_Data/RunC_5/plots/"

etas = [0.000, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.930,
        2.172, 2.322, 2.500, 2.650, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191]

etas_full = [-5.191, -3.839, -3.489, -3.139, -2.964, -2.853, -2.65, -2.5, -2.322, -2.172, -1.93, -1.653, -1.479, -1.305, -1.044, -0.783, -
             0.522, -0.261, 0, 0.261, 0.522, 0.783, 1.044, 1.305, 1.479, 1.653, 1.93, 2.172, 2.322, 2.5, 2.65, 2.853, 2.964, 3.139, 3.489, 3.839, 5.191]

etas_str = ["00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193",
            "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"]

etas_str_full = ["-5191", "-3839", "-3489", "-3139", "-2964", "-2853", "-265", "-25", "-2322", "-2172", "-193", "-1653", "-1479", "-1305", "-1044", "-0783", "-0522",
                 "-0261", "00", "0261", "0522", "0783", "1044", "1305", "1479", "1653", "193", "2172", "2322", "25", "2650", "2853", "2964", "3139", "3489", "3839", "5191"]

# pt_str_full = ["40", "62", "85", "153", "222", "290",
#                "359", "443", "486", "536", "1000", "2000"]

pt_str_full = ["51",
               "73",
               "85",
               "97",
               "179",
               "307",
               "370",
               "434",
               "520",
               "648", "1000", "2000", ]

pt_str = pt_str_full[1:]

trigger_val = ["40", "60", "80", "140",
               "200", "260", "320", "400",  # "450",
               "500"]

trigger_val_HF = ["60", "80", "100", "160",
                  "220", "300"]


pic_bases = [  # "alpha_pythia8",
    # "A_NormDistribution_DATA_pythia8",
    # "A_NormDistribution_MC_pythia8",
    "Rel_Response_pythia8",
    # "B_NormDistribution_DATA_pythia8",
    # "B_NormDistribution_MC_pythia8",
    #"jet3_pt_pythia8",
    #"METoverPt_DATA_pythia8",
    #"MET_pythia8",
    "MPF_Response_pythia8",
    # "Pt_ave_pythia8",

    # "probejet_phi_DATA_pythia8",
    # "probejet_phi_MC_pythia8",

    # "probejet_chEmEF_DATA_pythia8",
    # # "probejet_chEmEF_MC_pythia8",
    # "probejet_chHadEF_DATA_pythia8",
    # # "probejet_chHadEF_MC_pythia8",
    # "probejet_muonEF_DATA_pythia8",
    # # "probejet_muonEF_MC_pythia8",
    # "probejet_neutEmEF_DATA_pythia8",
    # # "probejet_neutEmEF_MC_pythia8",
    # "probejet_neutHadEF_DATA_pythia8",
    # # "probejet_neutHadEF_MC_pythia8",
    # "probejet_photonEF_DATA_pythia8",
    # # "probejet_photonEF_MC_pythia8",
]

pic_bases_MC = [  # "alpha_pythia8",
    # "A_NormDistribution_DATA_pythia8",
    "A_NormDistribution_MC_pythia8",
    # "Rel_Response_pythia8",
    # "B_NormDistribution_DATA_pythia8",
    "B_NormDistribution_MC_pythia8",
    #"jet3_pt_pythia8",
    #"METoverPt_DATA_pythia8",
    #"MET_pythia8",
    # "MPF_Response_pythia8",
    # "Pt_ave_pythia8",

    # "probejet_phi_DATA_pythia8",
    # "probejet_phi_MC_pythia8",

    # # "probejet_chEmEF_DATA_pythia8",
    # "probejet_chEmEF_MC_pythia8",
    # # "probejet_chHadEF_DATA_pythia8",
    # "probejet_chHadEF_MC_pythia8",
    # # "probejet_muonEF_DATA_pythia8",
    # "probejet_muonEF_MC_pythia8",
    # # "probejet_neutEmEF_DATA_pythia8",
    # "probejet_neutEmEF_MC_pythia8",
    # # "probejet_neutHadEF_DATA_pythia8",
    # "probejet_neutHadEF_MC_pythia8",
    # # "probejet_photonEF_DATA_pythia8",
    # "probejet_photonEF_MC_pythia8",
]


single_pics = [
    "kFSR_AK4PFchs_pythia8",
    "L2Res_kFSRfit_AK4PFchs_pythia8",
    "L2Res_kFSRval_AK4PFchs_pythia8",
    "L2Res_logpt_DiJet_kFSRfit_AK4PFchs_pythia8",
    # "L2Res_logpt_DiJet_kFSRfit_diffAK4PFchs_pythia8",
    "L2Res_logpt_MPF_DiJet_kFSRfitAK4PFchs_pythia8",
    "L2Res_logpt_MPF_kFSRfit_AK4PFchs_pythia8",
    # "L2Res_logpt_MPF_kFSRfit_diffAK4PFchs_pythia8",
]

first_pics = [
    "kFSR_AK4PFchs_pythia8",
    "L2Res_kFSRfit_AK4PFchs_pythia8",
    "L2Res_kFSRval_AK4PFchs_pythia8",
    "L2Res_logpt_MPF_kFSRfit_AK4PFchs_pythia8",
    "L2Res_logpt_DiJet_kFSRfit_AK4PFchs_pythia8",
    "L2Res_logpt_MPF_DiJet_kFSRfitAK4PFchs_pythia8",
]

document_base = """
\documentclass[t,compress]{beamer}
\mode<presentation>{
%\usetheme{CambridgeUS}
\usetheme{Madrid}
\usecolortheme{beaver}
}

\usepackage{graphicx} % Allows including images
\usepackage[percent]{overpic}
\usepackage{hyperref}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{upgreek}
%\usepackage{adjustbox}
\usepackage{setspace}

\usepackage{amsmath,amssymb}
\usepackage{color}

\usepackage{xcolor}
\usepackage{pict2e}

\usepackage{overpic}

\\newcommand{\\backupbegin}{
   \\newcounter{finalframe}
   \setcounter{finalframe}{\\value{framenumber}}
}
\\newcommand{\\backupend}{
   \setcounter{framenumber}{\\value{finalframe}}
}



\\title[L2 residual corrections]{2017 L2 residual corrections derived on di-jet events: State in September and Controlplots}
\\author[Christoph Garbers]{Anastasia Karavdina, Arne Reimers, Jens Multhaup, \underline{Christoph Garbers}}
\institute[UHH]{University of Hamburg}

\date{\\today}

\\begin{document}

\\begin{frame}
 \\titlepage
\end{frame}


"""

document_end = '\\backupend\n\end{document}'

frames_list = []

################ main frames #########################

frames_list.append(
    """
    \\begin{frame}{Differences to 2016 data}
    \\begin{itemize}
    \item For RunB and RunC no HLT Di-Jet Trigger available, Single-Jet Trigger (HLT PFJet)  followed by a Di-Jet selectin are used. \n\\newline\nUsed Data: JetH PromptReco RunB and RunC
    \item Compare to 2016 MC until 2017 is ready. \n\\newline\nUsed MC: Flat QCD pythia8
    \item Used JEC version: Fall17 17Fall2017 V4
    \end{itemize}
    \end{frame}

    """
)


################## phi vs eta ################

frames_list.append("\\begin{frame}{occupancy all jets Run D}\n")
frames_list[-1] += "\t\includegraphics[width=\\textwidth]{"
frames_list[-1] += img_path_base_D
frames_list[-1] += "ControlPlots_Selection_pythia8_PhiVsEtaAll.pdf}\n"
frames_list[-1] += "\end{frame}\n\n"

frames_list.append(
    "\\begin{frame}{occupancy probe jets Run D}\n")
frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_D
frames_list[-1] += "ControlPlots_Selection_pythia8_PhiVsEtaProbe.pdf}\n"
frames_list[-1] += "\end{frame}\n\n"


frames_list.append(
    "\\begin{frame}{occupancy all jets Run E}\n")
frames_list[-1] += "\t\includegraphics[width=0.8\\textwidth]{" + \
    img_path_base_E
frames_list[-1] += "ControlPlots_Selection_pythia8_PhiVsEtaAll.pdf}\n"
frames_list[-1] += "\end{frame}\n\n"

frames_list.append(
    "\\begin{frame}{occupancy probe jets Run E}\n")
frames_list[-1] += "\t\includegraphics[width=0.8\\textwidth]{" + \
    img_path_base_E
frames_list[-1] += "ControlPlots_Selection_pythia8_PhiVsEtaProbe.pdf}\n"
frames_list[-1] += "\end{frame}\n\n"

################## NPV vs eta ################


frames_list.append(
    "\\begin{frame}{NPV to eta Run D}\n")
for j, pt in enumerate(pt_str[:-1]):
    if j and not j % 5:
        frames_list[-1] += "\\newline\n\n"
    frames_list[-1] += "\t\includegraphics[width=0.19\\textwidth]{" + img_path_base_D + \
        "nvert_eta_DATA_pythia8_pt_" + \
        pt_str[j] + "_" + pt_str[j + 1] + ".pdf}\n"
frames_list[-1] += "\end{frame}\n\n"

frames_list.append(
    "\\begin{frame}{NPV to eta Run E}\n")
for j, pt in enumerate(pt_str[:-1]):
    if j and not j % 5:
        frames_list[-1] += "\\newline\n\n"
    frames_list[-1] += "\t\includegraphics[width=0.19\\textwidth]{" + img_path_base_E + \
        "nvert_eta_DATA_pythia8_pt_" + \
        pt_str[j] + "_" + pt_str[j + 1] + ".pdf}\n"
frames_list[-1] += "\end{frame}\n\n"


### Run  E #####
################ kFSR #########################


frames_list.append(
    "\\begin{frame}{L2Res with 2017 Run D}\n")
for i, spic in enumerate(first_pics):
    if i and not i % 3:
        frames_list[-1] += "\\newline\n\n"
    frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
    frames_list[-1] += "\t\\tiny  " +\
        spic.replace("_", " ").replace("pythia8", " ").replace(
            "AK4PFchs", " ").replace("logpt", " ") + "\n \t\\newline\n"
    frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_baseMain_D +\
        spic + ".pdf}\n"
    frames_list[-1] += "\end{minipage}\n"
frames_list[-1] += "\end{frame}\n\n"

frames_list.append(
    "\\begin{frame}{L2Res with 2017 Run E}\n")
for i, spic in enumerate(first_pics):
    if i and not i % 3:
        frames_list[-1] += "\\newline\n\n"
    frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
    frames_list[-1] += "\t\\tiny  " +\
        spic.replace("_", " ").replace("pythia8", " ").replace(
            "AK4PFchs", " ").replace("logpt", " ") + "\n \t\\newline\n"
    frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_baseMain_E +\
        spic + ".pdf}\n"
    frames_list[-1] += "\end{minipage}\n"
frames_list[-1] += "\end{frame}\n\n"


# frames_list.append(
#     "\\begin{frame}{L2Res with 2017 Run B}\n")
# for i, spic in enumerate(first_pics):
#     if i and not i % 3:
#         frames_list[-1] += "\\newline\n\n"
#     frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
#     frames_list[-1] += "\t\\tiny  " + \
#         spic.replace("_", " ").replace("pythia8", " ").replace(
#             "AK4PFchs", " ").replace("logpt", " ") + "\n \t\\newline\n"
#     frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_baseMain_B + \
#         spic + ".pdf}\n"
#     frames_list[-1] += "\end{minipage}\n"
# frames_list[-1] += "\end{frame}\n\n"

# frames_list.append(
#     "\\begin{frame}{L2Res with 2017 Run C}\n")
# for i, spic in enumerate(first_pics):
#     if i and not i % 3:
#         frames_list[-1] += "\\newline\n\n"
#     frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
#     frames_list[-1] += "\t\\tiny  " + \
#         spic.replace("_", " ").replace("pythia8", " ").replace(
#             "AK4PFchs", " ").replace("logpt", " ") + "\n \t\\newline\n"
#     frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_baseMain_C + \
#         spic + ".pdf}\n"
#     frames_list[-1] += "\end{minipage}\n"
# frames_list[-1] += "\end{frame}\n\n"


################ thresholds #########################

frames_list.append(
    "\\begin{frame}{PF DiJet Trigger Thresholds 0.95}\n")
frames_list[-1] += """
\\begin{tabular}{|c|c|c|}
  \hline
  trigger path & run DE & pre legacy 16 \\\\\hline
  HLT\_DiPFJet40\_v*             & 73 & 51 \\\\ \hline
  HLT\_DiPFJet60\_v*             & 85 & 73 \\\\ \hline
  HLT\_DiPFJet80\_v*             & 97 & 95 \\\\ \hline
  HLT\_DiPFJet140\_v*            & 179 & 129\\\\ \hline
  HLT\_DiPFJet200\_v*            & 307 & 163\\\\ \hline
  HLT\_DiPFJet260\_v*            & 370 & 230\\\\ \hline
  HLT\_DiPFJet320\_v*            & 434 & 299\\\\ \hline
  HLT\_DiPFJet400\_v*            & 520 & 365\\\\ \hline
  HLT\_DiPFJet500\_v*            & 649 & 566\\\\ \hline
  HLT\_DiPFJetAve60\_HFJEC\_v*  & 73  & \\\\ \hline
  HLT\_DiPFJetAve80\_HFJEC\_v*  & 93  & \\\\ \hline
  HLT\_DiPFJetAve100\_HFJEC\_v* & 113  & \\\\ \hline
  HLT\_DiPFJetAve160\_HFJEC\_v* & 176  & \\\\ \hline
  HLT\_DiPFJetAve220\_HFJEC\_v* & 239  & \\\\ \hline
  HLT\_DiPFJetAve300\_HFJEC\_v* & 318  & \\\\ \hline

\end{tabular}
 thresholds were derived without the jet matching for Run D and Run E seperatly
"""
frames_list[-1] += "\end{frame}\n\n"


# frames_list.append(
#     "\\begin{frame}{DiJet Trigger Thresholds center Run D}\n")
# frames_list[-1] += "\t\scriptsize shown is $\\varepsilon = \\frac{N( \\text{high trg})}{N(\\text{low trg})}$, dotted (full) line marks 0.9 (0.95) of plateau\n\t\\newline\nThe plateaus are not at one cause of different prescales of the different triggers.\n\\newline\n"

# for i, th in enumerate(trigger_val[1:]):
#     if i and not i % 4:
#         frames_list[-1] += "\t\\newline\n\n"
#     if not i == 4:
#         frames_list[-1] += "\\begin{minipage}{0.24\\textwidth}\n"
#         frames_list[-1] += "\t\\tiny HLT " + th + "\n \\newline\n"
#         frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" +\
#             img_path_base_th_D + "/HLT_DiPFJetAve"
#         frames_list[-1] += th
#         frames_list[-1] += ".pdf}\n"
#         frames_list[-1] += "\end{minipage}\n"

# frames_list[-1] += "\\begin{minipage}{0.24\\textwidth}\n"
# frames_list[-1] += "\t\\tiny extrapolation\n \\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" +\
#     img_path_base_th_D + "/extrapolateLowestTrigger.pdf}\n"
# frames_list[-1] += "\end{minipage}\n"

# frames_list[-1] += "\end{frame}\n\n"


# frames_list.append(
#     "\\begin{frame}{DiJet Trigger Thresholds fwd Run D}\n")
# frames_list[-1] += "\t\scriptsize shown is $\\varepsilon = \\frac{N( \\text{high trg})}{N(\\text{low trg})}$, dotted (full) line marks 0.9 (0.95) of plateau\n\t\\newline\nThe plateaus are not at one cause of different prescales of the different triggers.\n\\newline\n"

# for i, th in enumerate(trigger_val_HF[1:]):
#     if i and not i % 4:
#         frames_list[-1] += "\t\\newline\n\n"
#     if not i == 4:
#         frames_list[-1] += "\\begin{minipage}{0.24\\textwidth}\n"
#         frames_list[-1] += "\t\\tiny HLT " + th + "_HF\n \\newline\n"
#         frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" +\
#             img_path_base_th_D + "/HLT_DiPFJetAve"
#         frames_list[-1] += th
#         frames_list[-1] += ".pdf}\n"
#         frames_list[-1] += "\end{minipage}\n"

# frames_list[-1] += "\\begin{minipage}{0.24\\textwidth}\n"
# frames_list[-1] += "\t\\tiny extrapolation\n \\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" +\
#     img_path_base_th_D + "/extrapolateLowestTrigger_HF.pdf}\n"
# frames_list[-1] += "\end{minipage}\n"

# frames_list[-1] += "\end{frame}\n\n"


# frames_list.append(
#     "\\begin{frame}{DiJet Trigger Thresholds center Run E}\n")
# frames_list[-1] += "\t\scriptsize shown is $\\varepsilon = \\frac{N( \\text{high trg})}{N(\\text{low trg})}$, dotted (full) line marks 0.9 (0.95) of plateau\n\t\\newline\nThe plateaus are not at one cause of different prescales of the different triggers.\n\\newline\n"

# for i, th in enumerate(trigger_val[1:]):
#     if i and not i % 4:
#         frames_list[-1] += "\t\\newline\n\n"
#     if not i == 4:
#         frames_list[-1] += "\\begin{minipage}{0.24\\textwidth}\n"
#         frames_list[-1] += "\t\\tiny HLT " + th + "\n \\newline\n"
#         frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" +\
#             img_path_base_th_E + "/HLT_DiPFJetAve"
#         frames_list[-1] += th
#         frames_list[-1] += ".pdf}\n"
#         frames_list[-1] += "\end{minipage}\n"

# frames_list[-1] += "\\begin{minipage}{0.24\\textwidth}\n"
# frames_list[-1] += "\t\\tiny extrapolation\n \\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" +\
#     img_path_base_th_E + "/extrapolateLowestTrigger.pdf}\n"
# frames_list[-1] += "\end{minipage}\n"

# frames_list[-1] += "\end{frame}\n\n"

# frames_list.append(
#     "\\begin{frame}{DiJet Trigger Thresholds fwd Run E}\n")
# frames_list[-1] += "\t\scriptsize shown is $\\varepsilon = \\frac{N( \\text{high trg})}{N(\\text{low trg})}$, dotted (full) line marks 0.9 (0.95) of plateau\n\t\\newline\nThe plateaus are not at one cause of different prescales of the different triggers.\n\\newline\n"

# for i, th in enumerate(trigger_val_HF[1:]):
#     if i and not i % 4:
#         frames_list[-1] += "\t\\newline\n\n"
#     if not i == 4:
#         frames_list[-1] += "\\begin{minipage}{0.24\\textwidth}\n"
#         frames_list[-1] += "\t\\tiny HLT " + th + "_HF\n \\newline\n"
#         frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" +\
#             img_path_base_th_E + "/HLT_DiPFJetAve"
#         frames_list[-1] += th
#         frames_list[-1] += ".pdf}\n"
#         frames_list[-1] += "\end{minipage}\n"

# frames_list[-1] += "\\begin{minipage}{0.24\\textwidth}\n"
# frames_list[-1] += "\t\\tiny extrapolation\n \\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" +\
#     img_path_base_th_E + "/extrapolateLowestTrigger_HF.pdf}\n"
# frames_list[-1] += "\end{minipage}\n"

# frames_list[-1] += "\end{frame}\n\n"


########## Asym Examples ##########


# frames_list.append(
#     "\\begin{frame}{Not Understood: Asymmetry Distribution}\n")
# frames_list[-1] += "\t\scriptsize $A=\\frac{p_t^{tag}-p_t^{probe}}{p_t^{tag}+p_t^{probe}}$, jets are assigned randomly from the two leading jets. Only one trigger is used per pt-bin.\n\\newline\n"
# i = 0
# j = 2
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run B\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_B +\
#     "fullAsym/A_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# i = 3
# j = 8
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run B\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_B +\
#     "fullAsym/A_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# i = 9
# j = 8
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run B\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_B +\
#     "fullAsym/A_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# frames_list[-1] += "\\newline\n"

# frames_list[-1] += "\end{frame}\n\n"

# i = 0
# j = 2
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run C\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_C +\
#     "fullAsym/A_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# i = 3
# j = 8
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run C\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_C +\
#     "fullAsym/A_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# i = 9
# j = 8
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run C\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_C +\
#     "fullAsym/A_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\t\n"
# frames_list[-1] += "\end{frame}\n\n"


# frames_list.append(
#     "\\begin{frame}{Not Understood: Asymmetry Distribution, does not look pile-up dependend}\n")
# i = 0
# j = 2
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run B\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_B +\
#     "fullAsym/A_nvert_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# i = 3
# j = 8
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run B\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_B +\
#     "fullAsym/A_nvert_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# i = 9
# j = 8
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run B\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_B +\
#     "fullAsym/A_nvert_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# frames_list[-1] += "\\newline\n"

# i = 0
# j = 2
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run C\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_C +\
#     "fullAsym/A_nvert_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# i = 3
# j = 8
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run C\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_C +\
#     "fullAsym/A_nvert_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\n"
# i = 9
# j = 8
# frames_list[-1] += "\\begin{minipage}{0.3\\textwidth}\n"
# frames_list[-1] += "\t\scriptsize Run C\n\t\\newline\n"
# frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_C +\
#     "fullAsym/A_nvert_DATA_pythia8" + "_eta_" + etas_str[i] +\
#     "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
#     "_" + pt_str_full[j + 1] + ".pdf}\n"
# frames_list[-1] += "\end{minipage}\t\n"
# frames_list[-1] += "\end{frame}\n\n"


frames_list.append(
    "\\begin{frame}{Jet Properties after Selection for Run D}\n")
frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_D +\
    "ControlPlots_Selection_pythia8_jets.pdf} \n"
frames_list[-1] += "\end{frame}\n\n"

frames_list.append(
    "\\begin{frame}{Jet Properties after Selection for Run E}\n")
frames_list[-1] += "\t\includegraphics[width=\\textwidth]{" + img_path_base_E +\
    "ControlPlots_Selection_pythia8_jets.pdf} \n"
frames_list[-1] += "\end{frame}\n\n"

################ backup ##################################################

frames_list.append(
    """
    \\backupbegin
    \\begin{frame}
    Backup
    \end{frame}

    """
)


# eta loop D

for i, eta in enumerate(etas[:-1]):
    frames_list.append(
        "\\begin{frame}{Run D, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")

    for j, pic_base in enumerate(pic_bases[:4]):
        if j and not j % 2:
            frames_list[-1] += "\\newline\n\n"
        frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_D + "control/" +\
            pic_base + "_eta_" + etas_str[i] +\
            "_" + etas_str[i + 1] + ".pdf}\n"

    frames_list[-1] += "\end{frame}\n\n"
    # frames_list.append(
    #     "\\begin{frame}{Run D, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")

    # for j, pic_base in enumerate(pic_bases[4:10]):
    #     if j and not j % 3:
    #         frames_list[-1] += "\\newline\n\n"
    #         frames_list[-1] += "\\vspace{-0.65cm}\n"
    #     frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_D + "control/" +\
    #         pic_base + "_eta_" + etas_str[i] +\
    #         "_" + etas_str[i + 1] + ".pdf}\n"

    # frames_list[-1] += "\end{frame}\n\n"

    frames_list.append(
        "\\begin{frame}{Run D, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")
    frames_list[-1] += "\t Pt: \includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_D +\
        "kFSR_Pt" + "_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_D +\
        "kFSR_Pt_pythia8" + "_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\\newline\n\n"
    frames_list[-1] += "\t MPF: \includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_D +\
        "kFSR_MPF" + "_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_D +\
        "kFSR_MPF_pythia8" + "_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\end{frame}\n\n"

    frames_list.append(
        "\\begin{frame}{Run D Jet Pt, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_base_D +\
        "fullAsym/jet1_pt_DATA_pythia8_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_base_D +\
        "fullAsym/jet2_pt_DATA_pythia8_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_base_D +\
        "fullAsym/jet3_pt_DATA_pythia8_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\end{frame}\n\n"


# eta loop E

for i, eta in enumerate(etas[:-1]):
    frames_list.append(
        "\\begin{frame}{Run E, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")

    for j, pic_base in enumerate(pic_bases[:4]):
        if j and not j % 2:
            frames_list[-1] += "\\newline\n\n"
        frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_E + "control/" +\
            pic_base + "_eta_" + etas_str[i] +\
            "_" + etas_str[i + 1] + ".pdf}\n"

    frames_list[-1] += "\end{frame}\n\n"
    # frames_list.append(
    #     "\\begin{frame}{Run E, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")

    # for j, pic_base in enumerate(pic_bases[4:10]):
    #     if j and not j % 3:
    #         frames_list[-1] += "\\newline\n\n"
    #         frames_list[-1] += "\\vspace{-0.65cm}\n"
    #     frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_E + "control/" +\
    #         pic_base + "_eta_" + etas_str[i] +\
    #         "_" + etas_str[i + 1] + ".pdf}\n"

    # frames_list[-1] += "\end{frame}\n\n"

    frames_list.append(
        "\\begin{frame}{Run E, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")
    frames_list[-1] += "\t Pt: \includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_E +\
        "kFSR_Pt" + "_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_E +\
        "kFSR_Pt_pythia8" + "_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\\newline\n\n"
    frames_list[-1] += "\t MPF: \includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_E +\
        "kFSR_MPF" + "_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_baseMain_E +\
        "kFSR_MPF_pythia8" + "_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\end{frame}\n\n"

    frames_list.append(
        "\\begin{frame}{Run E Jet Pt, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_base_E +\
        "fullAsym/jet1_pt_DATA_pythia8_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_base_E +\
        "fullAsym/jet2_pt_DATA_pythia8_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\t\includegraphics[width=0.3\\textwidth]{" + img_path_base_E +\
        "fullAsym/jet3_pt_DATA_pythia8_eta_" + etas_str[i] +\
        "_" + etas_str[i + 1] + ".pdf}\n"
    frames_list[-1] += "\end{frame}\n\n"

# ############# eta loop MC  #########################

# for i, eta in enumerate(etas[:-1]):
#     frames_list.append(
#         "\\begin{frame}{MC, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")

#     for j, pic_base in enumerate(pic_bases_MC[:]):
#         if j and not j % 4:
#             frames_list[-1] += "\\newline\n\n"
#         frames_list[-1] += "\t\includegraphics[width=0.24\\textwidth]{" + img_path_baseMain_C + "control/" + \
#             pic_base + "_eta_" + etas_str[i] + \
#             "_" + etas_str[i + 1] + ".pdf}\n"

#     frames_list[-1] += "\end{frame}\n\n"


# #################### Asymmetry singled out #########################

for i, eta in enumerate(etas[:-1]):
    frames_list.append(
        "\\begin{frame}{Asymmetry Run D, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")
    for j, pt in enumerate(pt_str_full[:-1]):
        if j and not j % 5:
            frames_list[-1] += "\\newline\n\n"
        frames_list[-1] += "\t\includegraphics[width=0.19\\textwidth]{" + img_path_base_D +\
            "fullAsym/A_DATA_pythia8" + "_eta_" + etas_str[i] +\
            "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
            "_" + pt_str_full[j + 1] + ".pdf}\n"
    frames_list[-1] += "\end{frame}\n\n"


for i, eta in enumerate(etas[:-1]):
    frames_list.append(
        "\\begin{frame}{Asymmetry Run E, Eta Bin " + str(eta) + " - " + str(etas[i + 1]) + "}\n")
    for j, pt in enumerate(pt_str_full[:-1]):
        if j and not j % 5:
            frames_list[-1] += "\\newline\n\n"
        frames_list[-1] += "\t\includegraphics[width=0.19\\textwidth]{" + img_path_base_E +\
            "fullAsym/A_DATA_pythia8" + "_eta_" + etas_str[i] +\
            "_" + etas_str[i + 1] + "_pt_" + pt_str_full[j] +\
            "_" + pt_str_full[j + 1] + ".pdf}\n"
    frames_list[-1] += "\end{frame}\n\n"


#################### Put everything together #########################

with open(output_name, 'w') as file_:
    file_.write(document_base)
    for frame in frames_list:
        file_.write(frame)
    file_.write(document_end)

print "\ngenerated " + output_name + "\n"
