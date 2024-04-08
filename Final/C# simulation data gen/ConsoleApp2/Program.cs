using IsotopomerDynamics;
using RateConstant;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Data;
using System.Data.Common;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ConsoleApp2
{
    internal class Program
    {
        static void Main(string[] args)
        {
            #region other sim
            //#region generate sample mass isotopomers
            ////List<string> filecontent = new List<string>();
            ////filecontent.Add("M0,M1,M2,M3,M4,M5");
            ////for (double bwe = 0; bwe < 0.4; bwe= bwe+0.001)
            ////{
            ////    double[] fnatIso = new double[] { 38.37, 33.54, 17.92, 7.10, 2.27, 0.61 };
            ////    var res = NEH.generateTheroteicalIsotopomerValues((float)bwe, fnatIso.Select(x => (float)(x)).ToArray(), 24);
            ////    filecontent.Add(res[0] + "," + res[1] + "," + res[2] + "," + res[3] + "," + res[4] + "," + res[5]);
            ////}
            ////return;

            //#endregion

            //#region for neh simulation
            ///*
            //var filepath = @"C:\Workplace\Python\AnalysisForThePaper\NEH\second_round\prepare_Data\data_2.csv";
            //NEH neh_obj = new NEH();
            //var data = neh_obj.read_data(filepath);*/

            ///*
            //nEH.getExperimentalData((float)0.05, data);
            //nEH.getExperimentalData((float)0.04, data);
            //nEH.getExperimentalData((float)0.03, data);
            //nEH.getExperimentalData((float)0.02, data);
            //nEH.getExperimentalData((float)0.01, data);
            //*/

            ///*
            //#region neh simulation with different NEH NUMS
            //neh_obj.updateNEHValues(data);
            //data = data.Where(x => Math.Abs(x.M0 - x.I0) / x.M0 < 0.1).Take(3000).ToList();
            //neh_obj.getExperimentalData_fromExperimetntalI0((float)0.05, data);
            //neh_obj.computeNEHValuesBatch(data);
            //#endregion
            //*/

            //#region neh simulation with different NEH NUMS, error from io_asymp

            ////var filepath = @"C:\Workplace\Python\AnalysisForThePaper\NEH\second_round\prepare_Data\data_2.csv";
            ////NEH neh_obj = new NEH();
            ////var data = neh_obj.read_data(filepath);
            ////neh_obj.updateNEHValues(data);
            //////data = data.Where(x => Math.Abs(x.M0 - x.I0) / x.M0 < 0.1).Take(3000).ToList(); 
            ////data = neh_obj.getExperimentalData_fromExperimetntalI0_asymp((float)0.046, data);
            ////neh_obj.computeNEHValuesBatch(data, (float)0.046);

            ////var all_passed_peps = new List<PeptidesPassedNEHfilters>();
            ////foreach (var pep in data)
            ////{
            ////    PeptidesPassedNEHfilters peptidesPassedNEHfilters = new PeptidesPassedNEHfilters();
            ////    peptidesPassedNEHfilters.Peptide = pep.Peptide;
            ////    peptidesPassedNEHfilters.T_NEH = pep.NEH_sim;
            ////    peptidesPassedNEHfilters.MPE_NEH = pep.NEH_MPE;
            ////    peptidesPassedNEHfilters.MPE_RMSE = pep.MPE_rmse;
            ////    peptidesPassedNEHfilters.APE_NEH = pep.NEH_APE;
            ////    peptidesPassedNEHfilters.APE_RMSE = pep.APE_rmse;
            ////    all_passed_peps.Add(peptidesPassedNEHfilters);
            ////}
            ////List<AAsInfo> allAAs_info = NEH.getAllAAs();
            ////neh_obj.getNEHvalues(all_passed_peps, "MPE", allAAs_info);
            ////neh_obj.getNEHvalues(all_passed_peps, "APE", allAAs_info);

            //#endregion

            //#region two parameter test

            //var filepath = @"C:\Workplace\Python\AnalysisForThePaper\NEH\second_round\prepare_Data\data_3.csv";
            //NEH neh_obj = new NEH();
            //var data = neh_obj.read_data_withRIA(filepath);
            //neh_obj.updateNEHValues(data);
            //neh_obj.updateRIAVals_with_pxt_error(data);
            //////neh_obj.updateRIAVals_with_datapoint_noise(data);
            //neh_obj.compute_two_parm_rate_and_gof_neh(data);
            ////data = data.Where(x => x.two_param_rmse < 0.01 & x.two_param_rsquared > 0.97 & x.exp_times.Contains(31)).ToList();
            //data = data.Where(x => x.two_param_rmse < 0.1 & x.two_param_rsquared > 0.97 & x.exp_times.Contains(31)).ToList();
            ////data = data.Where(x => Math.Abs(x.RIAs.Last() - x.two_param_i0_asymp) / x.RIAs.Last() < 0.1).ToList();
            //data = data.Where(x => x.two_param_rate < Math.Log(2) & x.two_param_rate > Math.Log(2) / 31).ToList();
            ////data = data.Where(x => neh_obj.check_monotonically_decreasing_sequence(x.RIAs)).ToList();
            ////data = data.Where(x => Math.Abs(x.NEH_sim - x.NEH_two_param) > 10).ToList();

            //var all_passed_peps = new List<PeptidesPassedNEHfilters>();
            //foreach (var pep in data)
            //{
            //    PeptidesPassedNEHfilters peptidesPassedNEHfilters = new PeptidesPassedNEHfilters();
            //    peptidesPassedNEHfilters.Peptide = pep.Peptide;
            //    peptidesPassedNEHfilters.T_NEH = pep.NEH_sim;
            //    peptidesPassedNEHfilters.Asym_NEH = pep.NEH_two_param;
            //    peptidesPassedNEHfilters.Asym_RMSE = pep.two_param_rmse;
            //    all_passed_peps.Add(peptidesPassedNEHfilters);
            //}
            //List<AAsInfo> allAAs_info = NEH.getAllAAs();
            //neh_obj.getNEHvalues(all_passed_peps, "Asym", allAAs_info);

            //#endregion
            //#endregion


            #endregion

            #region rate constant test

            //Console.WriteLine("Test");

            //TDistributionCaller tdist = new TDistributionCaller();
            //Console.WriteLine((float)tdist.invcdf_T(0.025, 3 - 1));
            //Console.WriteLine((float)tdist.invcdf_T(0.975, 3 - 1));

            //var CurrentFolder = @"C:\Users\hmdebern.UTMB-USERS-M\Desktop\_abigail\test_lysate\working_3_1_1_2\\";
            //ProteinRateConstant ProtRate = new ProteinRateConstant(
            //    CurrentFolder,
            //    new float[] { 0, 0, 0, 0, 0, (float)(4 / 24.0), (float)4 / 24, (float)8 / 24, (float)8 / 24, (float)8 /24, 1, 1, 1, 2, 3 },
            //    new float[] { (float)0.04, (float)0.04, (float)0.04, (float)0.04, (float)0.04, 0, 0, 0, 0, 0,0,0,0,0,0},
            //    1, new List<float> { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
            //    (float)1.0, false, 2);
            //ProtRate.ReadFileFolder();

            ////var CurrentFolder = @"C:\Users\hmdebern.UTMB-USERS-M\Desktop\_abigail\test_lysate\lysate_output_3\\";
            //var CurrentFolder = @"C:\Users\hmdebern.UTMB-USERS-M\Desktop\_abigail\d2ome_output_2_1_25\lysate_output\\";
            //ProteinRateConstant ProtRate = new ProteinRateConstant(
            //    CurrentFolder,
            //    new float[] { 0, 0, 0, 0, 0, (float)(4 / 24.0), (float)4 / 24, (float)8 / 24, (float)8 / 24, (float)8 / 24, 1, 1, 1, 2, 3 },
            //    new float[] { (float)0.04, (float)0.04, (float)0.04, (float)0.04, (float)0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            //    1, new List<float> { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
            //    (float)1.0, false, 0);
            //ProtRate.ReadFileFolder();

            //var CurrentFolder = @"C:\Users\hmdebern.UTMB-USERS-M\Desktop\_abigail\test_lysate\lysate_output_3\\";
            //var CurrentFolder = @"C:\Users\hmdebern.UTMB-USERS-M\Desktop\_abigail\d2ome_output_2_1_25\lysate_output\\";
            //var CurrentFolder = @"C:\Users\hmdebern.UTMB-USERS-M\Desktop\_abigail\d2ome_output_2_1_25\hla_output\\";
            //ProteinRateConstant ProtRate = new ProteinRateConstant(
            //    CurrentFolder,
            //    new float[] { 0, 0, 0, 0, 0, (float)(4 / 24.0), (float)4 / 24, (float)8 / 24, (float)8 / 24, (float)8 / 24, 1, 1, 1, 2, 3 },
            //    new float[] { (float)0.04, (float)0.04, (float)0.04, (float)0.04, (float)0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            //    1, new List<float> { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
            //    (float)1.0, false, 0);
            //ProtRate.ReadFileFolder();


            ////DMSO
            //var CurrentFolder = @"C:\Users\hmdebern.UTMB-USERS-M\Desktop\_abigail\d2ome_output_4_2_24\DMSO\test2\\";
            //ProteinRateConstant ProtRate = new ProteinRateConstant(
            //    CurrentFolder,
            //    new float[] { 0, 0, 0, (float)(4 / 24.0), (float)4 / 24, (float)8 / 24, (float)8 / 24, 1, 1, 2, 2, 3 },
            //    new float[] { (float)0.04, (float)0.04, (float)0.04, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            //    1, new List<float> { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
            //    (float)1.0, false, 0);
            //ProtRate.ReadFileFolder();

            ////LEN
            //var CurrentFolder = @"C:\Users\hmdebern.UTMB-USERS-M\Desktop\_abigail\d2ome_output_4_2_24\LEN\test\\";
            //ProteinRateConstant ProtRate = new ProteinRateConstant(
            //    CurrentFolder,
            //    new float[] { 0, 0, 0, (float)(4 / 24.0), (float)8 / 24, 1, 2, 2, 2, 3 },
            //    new float[] { (float)0.04, (float)0.04, (float)0.04, 0, 0, 0, 0, 0, 0, 0 },
            //    1, new List<float> { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, },
            //    (float)1.0, false, 0);
            //ProtRate.ReadFileFolder();



            //var CurrentFolder = @"C:\Users\hmdebern.UTMB-USERS-M\Desktop\_abigail\test_lysate\working_1_1_1_2\\";
            //ProteinRateConstant ProtRate = new ProteinRateConstant(
            //    CurrentFolder,
            //    new float[] { 0,1,3 },
            //    new float[] { (float)0.04,0,0},
            //    1, new List<float> { 1, 1, 1},
            //    (float)1.0, false, 2);
            //ProtRate.ReadFileFolder();


            //ProteinRateConstant ^ ProtRate;
            //ProtRate = gcnew ProteinRateConstant(CurrentFolder, fTimeCourse, fBodyWaterEnrichment, Nparameter, abundance_scaling_factor,
            //    dTimeconversion, f_useTwo_mass_isotopomers, nPeptideConsistency);

            //ProtRate->ReadFileFolder();


            
            var CurrentFolder = @"C:\Workplace\executables_test\d2ome_4_1_2024\test_output\test\liverpool_liver\\";
            ProteinRateConstant ProtRate = new ProteinRateConstant(
                CurrentFolder,
                new float[] { 0, 1, 2, 3, 6, 7, 9, 13, 16, 21, 24, 31 },
                new float[] {(float)0,(float)0.03855948,(float)0.04479649,(float)0.04580533,(float)0.04599918,(float)0.04599987,
                    (float)0.046,(float)0.046,(float)0.046,(float)0.046,(float)0.046,(float)0.046},
                1, new List<float> { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
                (float)1.0, false, 0);
            ProtRate.ReadFileFolder();



            #endregion
        }
    }
}
