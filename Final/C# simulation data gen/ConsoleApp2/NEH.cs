using Accord.Math;
using Accord.Math.Optimization.Losses;
using Accord.Statistics.Models.Regression.Fitting;
using Accord.Statistics.Models.Regression.Linear;
using IsotopomerDynamics;
using LBFGS_Library_Call;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;


namespace ConsoleApp2
{
    public class NEH
    {

        public NEH()
        {

        }


        public List<DataRecord> read_data(string filepath)
        {

            List<DataRecord> dataRecords = new List<DataRecord>();

            {
                string[] lines = System.IO.File.ReadAllLines(filepath);
                lines = lines.Where(x => x.Length > 0).ToArray();
                for (int i = 1; i < lines.Length; i++)
                {
                    var fields = lines[i].Split(',');
                    // Assuming your CSV file has columns in the specified order
                    DataRecord record = new DataRecord
                    {
                        Protein = fields[0],
                        Peptide = fields[1],
                        Charge = (int)float.Parse(fields[2]), // Assuming Charge is an integer
                        NEH = int.Parse(fields[3]),    // Assuming NEH is an integer
                        M0 = float.Parse(fields[4]) / 100,     // Assuming M0 is an integer
                        M1 = float.Parse(fields[5]) / 100,     // Assuming M1 is an integer
                        M2 = float.Parse(fields[6]) / 100,     // Assuming M1 is an integer
                        M3 = float.Parse(fields[7]) / 100,     // Assuming M1 is an integer
                        M4 = float.Parse(fields[8]) / 100,     // Assuming M1 is an integer
                        M5 = float.Parse(fields[9]) / 100,     // Assuming M1 is an integer

                        I0 = float.Parse(fields[10]),     // Assuming M1 is an integer
                        I1 = float.Parse(fields[11]),     // Assuming M1 is an integer
                        I2 = float.Parse(fields[12]),     // Assuming M1 is an integer
                        I3 = float.Parse(fields[13]),     // Assuming M1 is an integer
                        I4 = float.Parse(fields[14]),     // Assuming M1 is an integer
                        I5 = float.Parse(fields[15]),     // Assuming M1 is an integer

                        I0_exp_t = float.Parse(fields[10 + 6]),     // Assuming M1 is an integer
                        I1_exp_t = float.Parse(fields[11 + 6]),     // Assuming M1 is an integer
                        I2_exp_t = float.Parse(fields[12 + 6]),     // Assuming M1 is an integer
                        I3_exp_t = float.Parse(fields[13 + 6]),     // Assuming M1 is an integer
                        I4_exp_t = float.Parse(fields[14 + 6]),     // Assuming M1 is an integer
                        I5_exp_t = float.Parse(fields[15 + 6]),     // Assuming M1 is an integer
                                                                    // Add other properties as needed
                    };

                    dataRecords.Add(record);
                }
            }

            return dataRecords;

        }



        public List<DataRecord> read_data_withRIA(string filepath)
        {

            List<DataRecord> dataRecords = new List<DataRecord>();

            {
                string[] lines = System.IO.File.ReadAllLines(filepath);
                lines = lines.Where(x => x.Length > 0).ToArray();
                for (int i = 1; i < lines.Length; i++)
                {
                    try
                    {
                        var fields = lines[i].Split(',');
                        // Assuming your CSV file has columns in the specified order
                        DataRecord record = new DataRecord
                        {
                            Protein = fields[0],
                            Peptide = fields[1],
                            Charge = (int)float.Parse(fields[2]), // Assuming Charge is an integer
                            NEH = int.Parse(fields[3]),    // Assuming NEH is an integer
                            M0 = float.Parse(fields[4]) / 100,     // Assuming M0 is an integer
                            M1 = float.Parse(fields[5]) / 100,     // Assuming M1 is an integer
                            M2 = float.Parse(fields[6]) / 100,     // Assuming M1 is an integer
                            M3 = float.Parse(fields[7]) / 100,     // Assuming M1 is an integer
                            M4 = float.Parse(fields[8]) / 100,     // Assuming M1 is an integer
                            M5 = float.Parse(fields[9]) / 100,     // Assuming M1 is an integer

                            I0 = float.Parse(fields[10]),     // Assuming M1 is an integer
                            I1 = float.Parse(fields[11]),     // Assuming M1 is an integer
                            I2 = float.Parse(fields[12]),     // Assuming M1 is an integer
                            I3 = float.Parse(fields[13]),     // Assuming M1 is an integer
                            I4 = float.Parse(fields[14]),     // Assuming M1 is an integer
                            I5 = float.Parse(fields[15]),     // Assuming M1 is an integer

                            I0_exp_t = float.Parse(fields[10 + 6]),     // Assuming M1 is an integer
                            I1_exp_t = float.Parse(fields[11 + 6]),     // Assuming M1 is an integer
                            I2_exp_t = float.Parse(fields[12 + 6]),     // Assuming M1 is an integer
                            I3_exp_t = float.Parse(fields[13 + 6]),     // Assuming M1 is an integer
                            I4_exp_t = float.Parse(fields[14 + 6]),     // Assuming M1 is an integer
                            I5_exp_t = float.Parse(fields[15 + 6])     // Assuming M1 is an integer
                                                                       // Add other properties as needed

                        };

                        List<int> exp_times = new List<int> { 0, 1, 2, 3, 6, 7, 9, 13, 16, 21, 24, 31 };
                        for (int index = 0; index < 12; index++)
                        {
                            float temp = 0;
                            float.TryParse(fields[22 + index], out temp);
                            if (temp != 0)
                            {
                                record.RIAs.Add(temp);
                                record.exp_times.Add(exp_times[index]);
                            }
                        }

                        //for (int index = 0; index < 12; index++)
                        //{
                        //    List<float> exp_ = new List<float>();
                        //    for (int index2 = 0; index2 < 6; index2++)
                        //    {
                        //        exp_.Add(float.Parse(fields[34 + (index * 6 + index2)]));
                        //    }
                        //    record.all_exp_spec.Add(exp_);
                        //}

                        dataRecords.Add(record);
                    }
                    catch (Exception e)
                    {
                        continue;
                    }
                }
            }

            return dataRecords;

        }

        public bool check_monotonically_decreasing_sequence(List<float> mylist)
        {
            for (int i = 0; i < mylist.Count - 1; i++)
            {
                if (mylist[i] < mylist[i + 1])
                    return false;
            }
            return true;
        }

        public void computeNEHValuesBatch(List<DataRecord> data, float pw = (float)0.05)
        {
            foreach (var record in data)
                compute_NEH(record, pw);
        }


        [DllImport("NNLS.dll", EntryPoint = "NonNegativeLeastSquares_3")]
        unsafe public static extern int NonNegativeLeastSquares_3(double[] a, int m, int n, double[] b, double[] x);

        public void compute_NEH(DataRecord peptide, float pw = (float)0.05)
        {

            #region APE

            List<double> current_Exp_isotope_distribution = new List<double>() { peptide.I0_exp_t, peptide.I1_exp_t, peptide.I2_exp_t,
                peptide.I3_exp_t, peptide.I4_exp_t, peptide.I5_exp_t, };
            current_Exp_isotope_distribution = normalize_list(current_Exp_isotope_distribution);



            var errors = new List<double>();
            var neh_list = new List<double>();
            var pxt_list = new List<double>();

            float[] fNatIsotopes2 = new float[6];
            float[] fLabIsotopes2 = new float[6];

            // compute the theoretical values
            for (double NEH_new2 = 1; NEH_new2 < 120; NEH_new2 = NEH_new2 + 0.1)
            {


                neh_list.Add(NEH_new2);
                pxt_list.Add(pw);

                fNatIsotopes2 = new float[6];
                fLabIsotopes2 = new float[6];
                fNatIsotopes2[0] = peptide.M0;
                fNatIsotopes2[1] = peptide.M1;
                fNatIsotopes2[2] = peptide.M2;
                fNatIsotopes2[3] = peptide.M3;
                fNatIsotopes2[4] = peptide.M4;
                fNatIsotopes2[5] = peptide.M5;
                MassIsotopomers MIDyn2 = new MassIsotopomers();


                MIDyn2.CalculateMIDynamics(fNatIsotopes2, fLabIsotopes2, pw, (float)NEH_new2);

                double[] theoretical_vals_at_t2 = { fLabIsotopes2[0], fLabIsotopes2[1], fLabIsotopes2[2], fLabIsotopes2[3],
                                                        fLabIsotopes2[4], fLabIsotopes2[5] };
                theoretical_vals_at_t2 = normalize_list(theoretical_vals_at_t2.ToList()).ToArray();

                var error = RMSE(current_Exp_isotope_distribution, theoretical_vals_at_t2.ToList());
                errors.Add(error);


            }


            var ape_neh = neh_list[errors.IndexOf(errors.Min())];
            var ape_neh_error = errors.Min();

            //if (ape_neh_error < 0.03)
            //    Console.WriteLine(peptide.NEH.ToString() + "\t" + ape_neh.ToString() + "\t" + Math.Abs(peptide.NEH - ape_neh).ToString());

            #endregion

            #region MPE method
            #region prepare "A" matrix


            double[][] A_vals = new double[6][];
            A_vals[0] = new double[6];
            A_vals[1] = new double[6];
            A_vals[2] = new double[6];
            A_vals[3] = new double[6];
            A_vals[4] = new double[6];
            A_vals[5] = new double[6];
            A_vals[0][0] = peptide.M0;
            A_vals[0][1] = peptide.M1;
            A_vals[0][2] = peptide.M2;
            A_vals[0][3] = peptide.M3;
            A_vals[0][4] = peptide.M4;
            A_vals[0][5] = peptide.M5;
            A_vals[1][0] = 0;
            A_vals[1][1] = A_vals[0][0];
            A_vals[1][2] = A_vals[0][1];
            A_vals[1][3] = A_vals[0][2];
            A_vals[1][4] = A_vals[0][3];
            A_vals[1][5] = A_vals[0][4];
            A_vals[2][0] = 0;
            A_vals[2][1] = 0;
            A_vals[2][2] = A_vals[0][0];
            A_vals[2][3] = A_vals[0][1];
            A_vals[2][4] = A_vals[0][2];
            A_vals[2][5] = A_vals[0][3];
            A_vals[3][0] = 0;
            A_vals[3][1] = 0;
            A_vals[3][2] = 0;
            A_vals[3][3] = A_vals[0][0];
            A_vals[3][4] = A_vals[0][1];
            A_vals[3][5] = A_vals[0][2];
            A_vals[4][0] = 0;
            A_vals[4][1] = 0;
            A_vals[4][2] = 0;
            A_vals[4][3] = 0;
            A_vals[4][4] = A_vals[0][0];
            A_vals[4][5] = A_vals[0][1];
            A_vals[5][0] = 0;
            A_vals[5][1] = 0;
            A_vals[5][2] = 0;
            A_vals[5][3] = 0;
            A_vals[5][4] = 0;
            A_vals[5][5] = A_vals[0][0];



            double[] a_array = new double[36];
            for (int q = 0; q < 6; q++)
            {
                for (int q1 = 0; q1 < 6; q1++)
                {
                    a_array[q * 6 + q1] = A_vals[q][q1];

                    //Console.WriteLine("" + (q * 6 + q1) + " =>" + q + " " + q1);
                }
            }


            double[] x_vals = new double[6];
            NonNegativeLeastSquares_3(a_array, 6, 6, current_Exp_isotope_distribution.ToArray(), x_vals);

            float px_t_new = (float)pw;
            double NEH_new = (((x_vals[2] / x_vals[1]) * 2 * ((1 - px_t_new) / px_t_new)) + 1);
            //double neh_x3 = ((x_vals[3] / x_vals[2]) * 3 * ((1 - px_t_new) / px_t_new)) + 2;

            #endregion


            #region compute neh error

            MassIsotopomers MIDyn = new MassIsotopomers();


            float[] fNatIsotopes = new float[6];
            float[] fLabIsotopes = new float[6];

            fNatIsotopes = new float[6];
            fNatIsotopes = new float[6];
            fNatIsotopes[0] = peptide.I0;
            fNatIsotopes[1] = peptide.I1;
            fNatIsotopes[2] = peptide.I2;
            fNatIsotopes[3] = peptide.I3;
            fNatIsotopes[4] = peptide.I4;
            fNatIsotopes[5] = peptide.I5;

            MIDyn.CalculateMIDynamics(fNatIsotopes, fLabIsotopes, px_t_new, (float)(NEH_new));//(int)Math.Round

            double[] theoretical_vals_at_t = { fLabIsotopes[0], fLabIsotopes[1], fLabIsotopes[2], fLabIsotopes[3], fLabIsotopes[4], fLabIsotopes[5] };
            theoretical_vals_at_t = normalize_list(theoretical_vals_at_t.ToList()).Select(x => (double)(x)).ToArray();
            var RMSE_error = RMSE(current_Exp_isotope_distribution.ToList(), theoretical_vals_at_t.ToList());

            //if (RMSE_error < 0.03)
            //    Console.WriteLine(peptide.NEH.ToString() + "\t" + NEH_new.ToString() + "\t" + Math.Abs(peptide.NEH - NEH_new).ToString());

            var mpe_neh = NEH_new;
            var mpe_neh_error = RMSE_error;
            #endregion

            #endregion

            peptide.NEH_APE = ape_neh;
            peptide.APE_rmse = (float)ape_neh_error;
            peptide.NEH_MPE = mpe_neh;// Math.Round(mpe_neh);
            peptide.MPE_rmse = (float)mpe_neh_error;

            //peptide.NEH_APE = peptide.NEH_sim;
            //peptide.APE_rmse = (float)0;
            //peptide.NEH_MPE = peptide.NEH_sim;
            //peptide.MPE_rmse = 0;

            if (peptide.Peptide == "AASDIAMTELPPTHPIR")
                Console.WriteLine("test");

        }


        public static double RMSE(List<double> selected_points, List<double> theoretical_points)
        {
            double rss = 0;
            //selected_points, theoretical_RIA
            for (int i = 0; i < selected_points.Count(); i++)
            {
                if (!double.IsNaN((double)(selected_points[i])))
                {
                    rss = rss + Math.Pow((double)(selected_points[i] - theoretical_points[i]), 2);
                }
            }
            var rmse = Math.Sqrt(rss / selected_points.Count());
            return rmse;
        }


        public void getExperimentalData_fromExperimetntalI0(float pw, List<DataRecord> data)
        {
            for (int index = 0; index < data.Count; index++)
            {
                #region compute theoretical isotopomer distribution 
                float[] _fNatIsotopes2 = new float[6];
                float[] _fLabIsotopes2 = new float[6];
                _fNatIsotopes2[0] = data[index].I0;
                _fNatIsotopes2[1] = data[index].I1;
                _fNatIsotopes2[2] = data[index].I2;
                _fNatIsotopes2[3] = data[index].I3;
                _fNatIsotopes2[4] = data[index].I4;
                _fNatIsotopes2[5] = data[index].I5;

                _fNatIsotopes2 = normalize_list(_fNatIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();

                MassIsotopomers _MIDyn = new MassIsotopomers();

                _MIDyn.CalculateMIDynamics(_fNatIsotopes2, _fLabIsotopes2, pw, (int)data[index].NEH_sim);

                _fLabIsotopes2 = normalize_list(_fLabIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();

                data[index].setExperimentalValues(_fLabIsotopes2);
                data[index].pxt = pw;
                #endregion
            }




        }


        public List<DataRecord> getExperimentalData_fromExperimetntalI0_asymp_original(float pw, List<DataRecord> data)
        {
            var res = new List<DataRecord>();

            for (int index = 0; index < data.Count; index++)
            {

                #region compute theoretical isotopomer distribution 
                float[] _fNatIsotopes2 = new float[6];
                float[] _fLabIsotopes2 = new float[6];
                _fNatIsotopes2[0] = data[index].I0;
                _fNatIsotopes2[1] = data[index].I1;
                _fNatIsotopes2[2] = data[index].I2;
                _fNatIsotopes2[3] = data[index].I3;
                _fNatIsotopes2[4] = data[index].I4;
                _fNatIsotopes2[5] = data[index].I5;
                _fNatIsotopes2 = normalize_list(_fNatIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();
                MassIsotopomers _MIDyn = new MassIsotopomers();
                _MIDyn.CalculateMIDynamics(_fNatIsotopes2, _fLabIsotopes2, pw, (float)(data[index].NEH));//(int)Math.Round
                _fLabIsotopes2 = normalize_list(_fLabIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();
                #endregion

                // experimental value at asymptote
                var peptide = data[index];
                List<double> current_Exp_isotope_distribution = new List<double>() { peptide.I0_exp_t, peptide.I1_exp_t, peptide.I2_exp_t,
                peptide.I3_exp_t, peptide.I4_exp_t, peptide.I5_exp_t, };
                current_Exp_isotope_distribution = normalize_list(current_Exp_isotope_distribution);

                var rmse_res = RMSE(
                    new List<double> { _fLabIsotopes2[0], _fLabIsotopes2[1], _fLabIsotopes2[2], _fLabIsotopes2[3], _fLabIsotopes2[4], _fLabIsotopes2[5] },
                    current_Exp_isotope_distribution
                    );
                if (rmse_res > 0.01) continue;


                // get error at the asymp.
                var errors = new List<double>();
                for (int i = 0; i < 6; i++)
                    errors.Add((_fLabIsotopes2[i] - current_Exp_isotope_distribution[i]) * 0.70);//* 1E-6


                #region recompute theoretical isotopomer distribution with sim. NEH and add noise
                _fNatIsotopes2 = new float[6];
                _fLabIsotopes2 = new float[6];
                _fNatIsotopes2[0] = data[index].M0;
                _fNatIsotopes2[1] = data[index].M1;
                _fNatIsotopes2[2] = data[index].M2;
                _fNatIsotopes2[3] = data[index].M3;
                _fNatIsotopes2[4] = data[index].M4;
                _fNatIsotopes2[5] = data[index].M5;
                _fNatIsotopes2 = normalize_list(_fNatIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();
                _MIDyn = new MassIsotopomers();
                _MIDyn.CalculateMIDynamics(_fNatIsotopes2, _fLabIsotopes2, pw, (float)(data[index].NEH_sim));//(int)Math.Round
                _fLabIsotopes2 = normalize_list(_fLabIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();


                var new_fLabIsotopes2 = new List<double>();
                for (int i = 0; i < 6; i++)
                    new_fLabIsotopes2.Add(_fLabIsotopes2[i] + errors[i]);
                new_fLabIsotopes2 = normalize_list(new_fLabIsotopes2);

                data[index].setExperimentalValues(new_fLabIsotopes2.Select(x => (float)x).ToArray());
                data[index].pxt = pw;
                #endregion
                res.Add(data[index]);
            }


            return res;


        }

        public List<DataRecord> getExperimentalData_fromExperimetntalI0_asymp(float pw, List<DataRecord> data)
        {
            var res = new List<DataRecord>();

            for (int index = 0; index < data.Count; index++)
            {

                // experimental value at asymptote
                //var peptide = data[index];


                var rmse_res = RMSE(
                    normalize_list(new List<double> { data[index].M0, data[index].M1, data[index].M2, data[index].M3, data[index].M4, data[index].M5 }),
                    normalize_list(new List<double> { data[index].I0, data[index].I1, data[index].I2, data[index].I3, data[index].I4, data[index].I5 })
                    );
                //Console.WriteLine(rmse_res);
                if (rmse_res > 0.1) continue;


                #region recompute theoretical isotopomer distribution with sim. NEH and add noise

                float[] _fNatIsotopes2 = new float[6];
                float[] _fLabIsotopes2 = new float[6];
                _fNatIsotopes2[0] = data[index].I0;
                _fNatIsotopes2[1] = data[index].I1;
                _fNatIsotopes2[2] = data[index].I2;
                _fNatIsotopes2[3] = data[index].I3;
                _fNatIsotopes2[4] = data[index].I4;
                _fNatIsotopes2[5] = data[index].I5;
                _fNatIsotopes2 = normalize_list(_fNatIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();
                var _MIDyn = new MassIsotopomers();
                _MIDyn.CalculateMIDynamics(_fNatIsotopes2, _fLabIsotopes2, pw, (float)(data[index].NEH_sim));//(int)Math.Round
                _fLabIsotopes2 = normalize_list(_fLabIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();


                data[index].setExperimentalValues(_fLabIsotopes2.Select(x => (float)x).ToArray());
                data[index].pxt = pw;
                #endregion
                res.Add(data[index]);
            }


            return res;


        }


        public static List<double> normalize_list(List<double> FirstArr)
        {
            //normalize first array
            var sum_FirstArray = FirstArr.Sum();
            var n_FirsArr = new List<double>();

            for (int i = 0; i < FirstArr.Count; i++)
            {
                n_FirsArr.Add(FirstArr[i] / sum_FirstArray);
            }
            return n_FirsArr;
        }

        public void getExperimentalData(float pw, List<DataRecord> data)
        {
            for (int index = 0; index < data.Count; index++)
            {
                #region compute theoretical isotopomer distribution 
                float[] _fNatIsotopes2 = new float[6];
                float[] _fLabIsotopes2 = new float[6];
                _fNatIsotopes2[0] = data[index].M0;// / 100;
                _fNatIsotopes2[1] = data[index].M1;// / 100;
                _fNatIsotopes2[2] = data[index].M2;// / 100;
                _fNatIsotopes2[3] = data[index].M3;// / 100;
                _fNatIsotopes2[4] = data[index].M4;// / 100;
                _fNatIsotopes2[5] = data[index].M5;// / 100;

                MassIsotopomers _MIDyn = new MassIsotopomers();

                _MIDyn.CalculateMIDynamics(_fNatIsotopes2, _fLabIsotopes2, pw, (int)data[index].NEH);

                data[index].setExperimentalValues(_fLabIsotopes2);
                data[index].pxt = pw;
                if (_fLabIsotopes2.Sum() < 0.6)
                    Console.WriteLine("=====>" + _fLabIsotopes2.Sum().ToString());
                #endregion
            }


            WriteCsvFile(pw.ToString().Replace('.', '_') + ".csv", data);



        }

        public void WriteCsvFile(string filePath, List<DataRecord> dataRecords)
        {

            TextWriter writer = new StreamWriter(filePath, false);

            // Write the header line
            writer.WriteLine("Protein,Peptide,Charge,NEH,M0,M1,M2,M3,M4,M5,I0,I1,I2,I3,I4,I5,pxt");
            foreach (var record in dataRecords)
            {
                string line = string.Format("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9},{10},{11},{12},{13},{14},{15},{16}",
                    record.Protein, record.Peptide, record.Charge, record.NEH, record.M0, record.M1,
                    record.M2, record.M3, record.M4, record.M5, record.I0_exp_t, record.I1_exp_t, record.I2_exp_t,
                    record.I3_exp_t, record.I4_exp_t, record.I5_exp_t, record.pxt);
                writer.WriteLine(line);
            }
            writer.Close();



            //using (StreamWriter writer = new StreamWriter(filePath, false))
            //{
            //    // Write the header line
            //    writer.WriteLine("Protein,Peptide,Charge,NEH,M0,M1,,M2,M3,M4,M5,I0,I1,I2,I3,I4,I5");

            //    // Write data for each DataRecord
            //    foreach (var record in dataRecords)
            //    {
            //        var line = $"{record.Protein},{record.Peptide},{record.Charge},{record.NEH},{record.M0},{record.M1},{record.M2},{record.M3},{record.M4},{record.M5},{record.I1},{record.I2},{record.I3},{record.I4},{record.I5}";
            //        writer.WriteLine(line);
            //    }
            //}
        }

        public void updateNEHValues(List<DataRecord> data)
        {
            foreach (var record in data)
                record.NEH_sim = getSimulationNEH(record.Peptide);
        }
        public static float getSimulationNEH(string seq)
        {

            //var aa_neh = new Dictionary<string, float>() { {"a", 4},{"c", 2},{"d", 2},{"e", 4},{"f", 1},{"g", 2},{"h", 3},
            //    {"i", 1},{"k", 1},{"l", 1},{"m", 1},{"n", 2},{"p", 3},{"q", 4},{"r", 4},{"s", 3},{"t", 0},{"v", 1},{"w", 0},{"y", 1}};
            //float neh = 0;
            //foreach (var aa in seq.ToLower())
            //{
            //    if (aa.Equals('k'))
            //        neh += (float)4.5;
            //    else if (aa.Equals('r'))
            //        neh += (float)2.3;
            //    else if (aa.Equals('a'))
            //        neh += (float)4.3;
            //    else if (aa.Equals('c'))
            //        neh += (float)2.3;
            //    else if (aa.Equals('d'))
            //        neh += (float)2.3;
            //    else if (aa.Equals('f'))
            //        neh += (float)3.3;
            //    else
            //        neh += (float)(1.25 * aa_neh[aa.ToString()]);
            //}
            //return neh;

            //var aa_neh = new Dictionary<string, double>() { {"a", 4},{"c", 1.62},{"d", 1.89},{"e", 3.95},{"f", 0.32},{"g", 2.06},{"h", 2.88},
            //    {"i", 1},{"k", 0.54},{"l", 0.6},{"m", 1.12},{"n", 1.89},{"p", 2.59},{"q", 3.95},{"r", 3.43},{"s", 2.61},{"t", 0.2},{"v", 0.56},{"w", 0.08},{"y", 0.42}};

            var aa_neh = new Dictionary<string, double>() { {"a", 2},{"c", 3},{"d", 3},{"e", 2},{"f", 1},{"g", 4},{"h", 1},
                {"i", 2},{"k", 2},{"l", 3},{"m", 0},{"n", 1},{"p", 0},{"q", 2},{"r", 1},{"s", 1},{"t",4},{"v", 0},{"w",3},{"y", 1}};


            double neh = 0;
            foreach (var aa in seq.ToLower())
                neh += (aa_neh[aa.ToString()]);
            return (float)neh;

        }


        public void updateRIAVals_with_pxt_error(List<DataRecord> data)
        {
            double ph = 1.5574E-4;
            foreach (var record in data)
            {
                for (int index = 0; index < record.RIAs.Count(); index++)
                {
                    if (record.RIAs[index] != 0)
                    {
                        var new_pxt = (1 - ph) * (1 - Math.Pow(record.RIAs[index] / record.M0, 1.0 / record.NEH));

                        float newr_ria = (float)(record.M0 * Math.Pow(1 - (new_pxt / (1 - ph)), record.NEH_sim));

                        record.RIAs[index] = newr_ria;
                        record.pxt = (float)new_pxt;
                        //if (index == 11)
                        //    Console.WriteLine(record.exp_times[index] + "," + new_pxt);
                    }
                }
            }

        }

        public void updateRIAVals_with_datapoint_noise(List<DataRecord> data)
        {
            double ph = 1.5574E-4;
            var pw = (float)0.046;
            foreach (var record in data)
            {
                /*
                for (int index = 0; index < record.RIAs.Count(); index++)
                {
                    if (record.RIAs[index] != 0)
                    {
                        var new_pxt = (1 - ph) * (1 - Math.Pow(record.RIAs[index] / record.M0, 1.0 / record.NEH));

                        float newr_ria = (float)(record.M0 * Math.Pow(1 - (new_pxt / (1 - ph)), record.NEH_sim));

                        record.RIAs[index] = newr_ria;
                        record.pxt = (float)new_pxt;
                        //if (index == 11)
                        //    Console.WriteLine(record.exp_times[index] + "," + new_pxt);
                    }
                }*/

                for (int timepointindex = 0; timepointindex < 12; timepointindex++)
                {
                    var current_Exp_isotope_distribution = normalize_list(record.all_exp_spec[timepointindex].Select(x => (double)(x)).ToList());
                    #region get theoretical spec

                    float[] _fNatIsotopes2 = new float[6];
                    float[] _fLabIsotopes2 = new float[6];
                    _fNatIsotopes2[0] = record.M0;
                    _fNatIsotopes2[1] = record.M1;
                    _fNatIsotopes2[2] = record.M2;
                    _fNatIsotopes2[3] = record.M3;
                    _fNatIsotopes2[4] = record.M4;
                    _fNatIsotopes2[5] = record.M5;
                    _fNatIsotopes2 = normalize_list(_fNatIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();
                    MassIsotopomers _MIDyn = new MassIsotopomers();
                    _MIDyn.CalculateMIDynamics(_fNatIsotopes2, _fLabIsotopes2, pw, (int)record.NEH);
                    _fLabIsotopes2 = normalize_list(_fLabIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();


                    //var rmse_res = RMSE(
                    //        new List<double> { _fLabIsotopes2[0], _fLabIsotopes2[1], _fLabIsotopes2[2], _fLabIsotopes2[3], _fLabIsotopes2[4], _fLabIsotopes2[5] },
                    //        current_Exp_isotope_distribution.Select(x => (double)x).ToList()
                    //        ); ;
                    //if (rmse_res > 0.01) continue;



                    #endregion


                    // get error 
                    var errors = new List<double>();
                    for (int i = 0; i < 6; i++)
                        errors.Add(_fLabIsotopes2[i] - current_Exp_isotope_distribution[i]);

                    #region recompute theoretical isotopomer distribution with sim. NEH and add noise
                    _fNatIsotopes2 = new float[6];
                    _fLabIsotopes2 = new float[6];
                    _fNatIsotopes2[0] = record.M0;
                    _fNatIsotopes2[1] = record.M1;
                    _fNatIsotopes2[2] = record.M2;
                    _fNatIsotopes2[3] = record.M3;
                    _fNatIsotopes2[4] = record.M4;
                    _fNatIsotopes2[5] = record.M5;
                    _fNatIsotopes2 = normalize_list(_fNatIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();
                    _MIDyn = new MassIsotopomers();
                    _MIDyn.CalculateMIDynamics(_fNatIsotopes2, _fLabIsotopes2, pw, (int)record.NEH_sim);
                    _fLabIsotopes2 = normalize_list(_fLabIsotopes2.Select(x => (double)(x)).ToList()).Select(x => (float)(x)).ToArray();


                    var new_fLabIsotopes2 = new List<double>();
                    for (int i = 0; i < 6; i++)
                        new_fLabIsotopes2.Add(_fLabIsotopes2[i] + errors[i]);

                    record.all_exp_spec[timepointindex] = normalize_list(new_fLabIsotopes2).Select(x => (float)x).ToList();
                    record.RIAs[timepointindex] = record.all_exp_spec[timepointindex][0];
                    #endregion
                }

            }

        }


        public void getNEHvalues(List<PeptidesPassedNEHfilters> passedPeptides, string Method, List<AAsInfo> allAAs_info)
        {
            float rmse_th = float.PositiveInfinity;// (float)0.01;
            //float rmse_th = float.MaxValue; //(float)numericUpDown_neh_min_specRMSE.Value; float rmse_th = (float)numericUpDown_neh_min_specRMSE.Value;
            List<PeptidesPassedNEHfilters> items = new List<PeptidesPassedNEHfilters>();
            List<double> neh_dif;
            int total_eq = 1000;
            //passedPeptides = passedPeptides.Where(x => x.Peptide.Contains('R')).ToList();
            //passedPeptides = passedPeptides.Where(x => x.Peptide.Length > 15 & x.Peptide.Length < 45).ToList();
            switch (Method)
            {
                case "Asym":
                    items = passedPeptides.Where(x => x.Asym_RMSE <= rmse_th).ToList();
                    //items = items.OrderBy(x => x.Asym_RMSE).Take(total_eq).ToList();
                    neh_dif = items.Select(x => (double)(100.0 * Math.Abs(x.T_NEH - x.Asym_NEH) / x.T_NEH)).ToList();

                    break;
                case "APE":
                    items = passedPeptides.Where(x => x.APE_RMSE <= rmse_th).ToList();
                    //items = items.OrderBy(x => x.APE_RMSE).Take(total_eq).ToList();
                    neh_dif = items.Select(x => (double)(100.0 * Math.Abs(x.T_NEH - x.APE_NEH) / x.T_NEH)).ToList();

                    break;
                case "MPE":
                    items = passedPeptides.Where(x => x.MPE_RMSE <= rmse_th).ToList();
                    //items = items.OrderBy(x => x.MPE_RMSE).Take(total_eq).ToList();
                    neh_dif = items.Select(x => (double)(100.0 * Math.Abs(x.T_NEH - x.MPE_NEH) / x.T_NEH)).ToList();

                    break;
                case "APE_MPE":
                    items = passedPeptides.Where(x => Math.Max(x.APE_NEH, Math.Max(x.MPE_NEH, x.Asym_NEH)) - Math.Min(x.APE_NEH, Math.Min(x.MPE_NEH, x.Asym_NEH)) < 2 &
                    x.APE_RMSE <= rmse_th & x.MPE_RMSE <= rmse_th & x.Asym_RMSE <= rmse_th).ToList();

                    neh_dif = items.Select(x => (double)(100.0 * Math.Abs(x.T_NEH - (x.MPE_NEH + x.Asym_NEH + x.APE_NEH) / 3) / x.T_NEH)).ToList();

                    break;
            }


            // this two variables are used to compute NNLS
            double[][] result_accord = new double[items.Count][];

            List<double> result = new List<double>(new double[items.Count * 20]);
            List<double> result_neh = new List<double>();

            // this two variables are used to compute CI
            Vector<double> true_neh = Vector<double>.Build.Dense(items.Count);
            Matrix<double> coef_matrix = Matrix<double>.Build.Random(items.Count, 20);

            for (int i = 0; i < items.Count; i++)
            {
                Dictionary<char, int> current_pepinfo = new Dictionary<char, int>();
                //foreach (AAsInfo aasInfo in allAAs_info)
                for (int aa_index = 0; aa_index < 20; aa_index++)
                {
                    current_pepinfo[allAAs_info[aa_index].AA_symbol.Trim()[0]] = 0;
                    coef_matrix[i, aa_index] = 0;
                }

                foreach (char c in items[i].Peptide.ToUpper())
                {
                    current_pepinfo[c] += 1;
                    coef_matrix[i, allAAs_info.IndexOf(allAAs_info.Where(x => x.AA_symbol.Trim()[0] == c).FirstOrDefault())] += 1.0;

                }


                List<double> temp_res = new List<double>();
                double[] accord_row = new double[20];
                for (int index_c = 0; index_c < 20; index_c++)
                {
                    AAsInfo aasInfo = allAAs_info[index_c];
                    result[i + items.Count * index_c] = current_pepinfo[aasInfo.AA_symbol.Trim()[0]];
                    accord_row[index_c] = current_pepinfo[aasInfo.AA_symbol.Trim()[0]];
                }
                result_accord[i] = accord_row;

                switch (Method)
                {
                    case "Asym":
                        result_neh.Add(items[i].Asym_NEH);
                        true_neh[i] = items[i].Asym_NEH;
                        break;
                    case "APE":
                        result_neh.Add(items[i].APE_NEH);
                        true_neh[i] = items[i].APE_NEH;
                        break;
                    case "MPE":
                        result_neh.Add(items[i].MPE_NEH);
                        true_neh[i] = items[i].MPE_NEH;
                        break;
                    case "APE_MPE":
                        result_neh.Add((items[i].MPE_NEH + items[i].APE_NEH + items[i].Asym_NEH) / 3);
                        true_neh[i] = ((items[i].MPE_NEH + items[i].APE_NEH + items[i].Asym_NEH) / 3);
                        break;
                }


            }

            double[] x_vals = new double[20];
            NonNegativeLeastSquares_3(result.ToArray(), items.Count, 20, result_neh.ToArray(), x_vals);

            #region accord statistics
            {
                // Declare training samples
                var inputs = result_accord;

                var outputs = result_neh.ToArray();

                // Create a NN LS learning algorithm
                var nnls = new NonNegativeLeastSquares()
                {
                    MaxIterations = 1000000
                };

                // Use the algorithm to learn a multiple linear regression
                MultipleLinearRegression regression = nnls.Learn(inputs, outputs);

                // None of the regression coefficients should be negative:
                double[] coefficients = regression.Weights; // should be

                // Check the quality of the regression:
                double[] prediction = regression.Transform(inputs);

                double error = new SquareLoss(expected: outputs)
                    .Loss(actual: prediction); // should be 0

                //x_vals = coefficients;

                //var temp_x_vals = (coef_matrix.Transpose() * coef_matrix).Inverse() * coef_matrix.Transpose() * true_neh;
                //for (int i = 0; i < 20; i++)
                //    x_vals[i] = temp_x_vals[i];

            }

            #endregion

            #region compute CI
            Vector<double> parms = Vector<double>.Build.Dense(20);
            Vector<double> parms_ci = Vector<double>.Build.Dense(20);
            for (int aa_index = 0; aa_index < 20; aa_index++)
                parms[aa_index] = x_vals[aa_index];

            var pred_neh = coef_matrix * parms;

            var resduals = ((true_neh - pred_neh) * (true_neh - pred_neh)) / (items.Count - 20);

            var tval = StudentT.InvCDF(location: 0.0, scale: 1.0, freedom: (items.Count - 20), p: 1.0 - 0.05 / 2.0);

            var cov = (coef_matrix.Transpose() * coef_matrix).Inverse();

            for (int aa_index = 0; aa_index < 20; aa_index++)
            {
                var vari = resduals * cov[aa_index, aa_index];
                parms_ci[aa_index] = tval * Math.Sqrt(vari);
            }

            #endregion

            for (int i = 0; i < x_vals.Length; i++)
            {
                switch (Method)
                {
                    case "Asym":
                        allAAs_info[i].AA_NEH_Asym = x_vals[i];
                        allAAs_info[i].AA_NEH_Asym_ci = parms_ci[i];
                        break;
                    case "APE":
                        allAAs_info[i].AA_NEH_APE = x_vals[i];
                        allAAs_info[i].AA_NEH_APE_ci = parms_ci[i];
                        break;
                    case "MPE":
                        allAAs_info[i].AA_NEH_MPE = x_vals[i];
                        allAAs_info[i].AA_NEH_MPE_ci = parms_ci[i];
                        break;
                    case "APE_MPE":
                        allAAs_info[i].AA_NEH_APE_MPE = x_vals[i];
                        allAAs_info[i].AA_NEH_APE_MPE_ci = parms_ci[i];
                        break;
                }

            }
        }


        public void getNEHvalues_original(List<PeptidesPassedNEHfilters> passedPeptides, string Method, List<AAsInfo> allAAs_info)
        {
            float rmse_th = (float)0.01;
            //float rmse_th = float.MaxValue; //(float)numericUpDown_neh_min_specRMSE.Value; float rmse_th = (float)numericUpDown_neh_min_specRMSE.Value;
            List<PeptidesPassedNEHfilters> items = new List<PeptidesPassedNEHfilters>();
            List<double> neh_dif;
            int total_eq = 1000;
            //passedPeptides = passedPeptides.Where(x => x.Peptide.Contains('R')).ToList();
            switch (Method)
            {
                case "Asym":
                    items = passedPeptides.Where(x => x.Asym_RMSE <= rmse_th).ToList();
                    //items = items.OrderBy(x => x.Asym_RMSE).Take(total_eq).ToList();
                    neh_dif = items.Select(x => (double)(100.0 * Math.Abs(x.T_NEH - x.Asym_NEH) / x.T_NEH)).ToList();

                    break;
                case "APE":
                    items = passedPeptides.Where(x => x.APE_RMSE <= rmse_th).ToList();
                    //items = items.OrderBy(x => x.APE_RMSE).Take(total_eq).ToList();
                    neh_dif = items.Select(x => (double)(100.0 * Math.Abs(x.T_NEH - x.APE_NEH) / x.T_NEH)).ToList();

                    break;
                case "MPE":
                    items = passedPeptides.Where(x => x.MPE_RMSE <= rmse_th).ToList();
                    //items = items.OrderBy(x => x.MPE_RMSE).Take(total_eq).ToList();
                    neh_dif = items.Select(x => (double)(100.0 * Math.Abs(x.T_NEH - x.MPE_NEH) / x.T_NEH)).ToList();

                    break;
                case "APE_MPE":
                    items = passedPeptides.Where(x => Math.Max(x.APE_NEH, Math.Max(x.MPE_NEH, x.Asym_NEH)) - Math.Min(x.APE_NEH, Math.Min(x.MPE_NEH, x.Asym_NEH)) < 2 &
                    x.APE_RMSE <= rmse_th & x.MPE_RMSE <= rmse_th & x.Asym_RMSE <= rmse_th).ToList();

                    neh_dif = items.Select(x => (double)(100.0 * Math.Abs(x.T_NEH - (x.MPE_NEH + x.Asym_NEH + x.APE_NEH) / 3) / x.T_NEH)).ToList();

                    break;
            }


            // this two variables are used to compute NNLS
            List<double> result = new List<double>(new double[items.Count * 20]);
            List<double> result_neh = new List<double>();

            // this two variables are used to compute CI
            Vector<double> true_neh = Vector<double>.Build.Dense(items.Count);
            Matrix<double> coef_matrix = Matrix<double>.Build.Random(items.Count, 20);

            for (int i = 0; i < items.Count; i++)
            {
                Dictionary<char, int> current_pepinfo = new Dictionary<char, int>();
                //foreach (AAsInfo aasInfo in allAAs_info)
                for (int aa_index = 0; aa_index < 20; aa_index++)
                {
                    current_pepinfo[allAAs_info[aa_index].AA_symbol.Trim()[0]] = 0;
                    coef_matrix[i, aa_index] = 0;
                }

                foreach (char c in items[i].Peptide.ToUpper())
                {
                    current_pepinfo[c] += 1;
                    coef_matrix[i, allAAs_info.IndexOf(allAAs_info.Where(x => x.AA_symbol.Trim()[0] == c).FirstOrDefault())] += 1;
                }


                List<double> temp_res = new List<double>();
                for (int index_c = 0; index_c < 20; index_c++)
                {
                    AAsInfo aasInfo = allAAs_info[index_c];
                    result[i + items.Count * index_c] = current_pepinfo[aasInfo.AA_symbol.Trim()[0]];
                }


                switch (Method)
                {
                    case "Asym":
                        result_neh.Add(items[i].Asym_NEH);
                        true_neh[i] = items[i].Asym_NEH;
                        break;
                    case "APE":
                        result_neh.Add(items[i].APE_NEH);
                        true_neh[i] = items[i].APE_NEH;
                        break;
                    case "MPE":
                        result_neh.Add(items[i].MPE_NEH);
                        true_neh[i] = items[i].MPE_NEH;
                        break;
                    case "APE_MPE":
                        result_neh.Add((items[i].MPE_NEH + items[i].APE_NEH + items[i].Asym_NEH) / 3);
                        true_neh[i] = ((items[i].MPE_NEH + items[i].APE_NEH + items[i].Asym_NEH) / 3);
                        break;
                }


            }

            double[] x_vals = new double[20];
            NonNegativeLeastSquares_3(result.ToArray(), items.Count, 20, result_neh.ToArray(), x_vals);

            #region compute CI
            Vector<double> parms = Vector<double>.Build.Dense(20);
            Vector<double> parms_ci = Vector<double>.Build.Dense(20);
            for (int aa_index = 0; aa_index < 20; aa_index++)
                parms[aa_index] = x_vals[aa_index];

            var pred_neh = coef_matrix * parms;

            var resduals = ((true_neh - pred_neh) * (true_neh - pred_neh)) / (items.Count - 20);

            var tval = StudentT.InvCDF(location: 0.0, scale: 1.0, freedom: (items.Count - 20), p: 1.0 - 0.05 / 2.0);

            var cov = (coef_matrix.Transpose() * coef_matrix).Inverse();

            for (int aa_index = 0; aa_index < 20; aa_index++)
            {
                var vari = resduals * cov[aa_index, aa_index];
                parms_ci[aa_index] = tval * Math.Sqrt(vari);
            }

            #endregion

            for (int i = 0; i < x_vals.Length; i++)
            {
                switch (Method)
                {
                    case "Asym":
                        allAAs_info[i].AA_NEH_Asym = x_vals[i];
                        allAAs_info[i].AA_NEH_Asym_ci = parms_ci[i];
                        break;
                    case "APE":
                        allAAs_info[i].AA_NEH_APE = x_vals[i];
                        allAAs_info[i].AA_NEH_APE_ci = parms_ci[i];
                        break;
                    case "MPE":
                        allAAs_info[i].AA_NEH_MPE = x_vals[i];
                        allAAs_info[i].AA_NEH_MPE_ci = parms_ci[i];
                        break;
                    case "APE_MPE":
                        allAAs_info[i].AA_NEH_APE_MPE = x_vals[i];
                        allAAs_info[i].AA_NEH_APE_MPE_ci = parms_ci[i];
                        break;
                }

            }
        }
        //compute_two_parm_rate_and_gof_neh

        #region AAs neh nums
        public static AAsInfo Alanine = new AAsInfo("Alanine (A)", "a", getSimulationNEH("a")); //4.0
        public static AAsInfo Cysteine = new AAsInfo("Cysteine (C)", "c", getSimulationNEH("c")); //1.62
        public static AAsInfo Aspartic_acid = new AAsInfo("Aspartic acid (D)", "d", getSimulationNEH("d")); //1.89);
        public static AAsInfo Glutamic_acid = new AAsInfo("Glutamic acid (E)", "e", getSimulationNEH("e")); // 3.95);
        public static AAsInfo Phenylalanine = new AAsInfo("Phenylalanine (F)", "f", getSimulationNEH("f")); //0.32);
        public static AAsInfo Glycine = new AAsInfo("Glycine (G)", "g", getSimulationNEH("g")); //2.06);
        public static AAsInfo Histidine = new AAsInfo("Histidine (H)", "h", getSimulationNEH("h")); //2.88);
        public static AAsInfo Isoleucine = new AAsInfo("Isoleucine (I)", "i", getSimulationNEH("i")); //1.0);
        public static AAsInfo lysine = new AAsInfo("lysine (K)", "k", getSimulationNEH("k")); //0.6);
        public static AAsInfo Leucine = new AAsInfo("Leucine (L)", "l", getSimulationNEH("l")); //0.54);
        public static AAsInfo Methionine = new AAsInfo("Methionine (M)", "m", getSimulationNEH("m")); //1.12);
        public static AAsInfo Asparagine = new AAsInfo("Asparagine (N)", "n", getSimulationNEH("n")); //1.89);
        public static AAsInfo Proline = new AAsInfo("Proline (P)", "p", getSimulationNEH("p")); //2.59);
        public static AAsInfo Glutamine = new AAsInfo("Glutamine (Q)", "q", getSimulationNEH("q")); //3.95);
        public static AAsInfo Arginine = new AAsInfo("Arginine (R)", "r", getSimulationNEH("r")); //3.43);
        public static AAsInfo Serine = new AAsInfo("Serine (S)", "s", getSimulationNEH("s")); //2.61);
        public static AAsInfo Threonine = new AAsInfo("Threonine (T)", "t", getSimulationNEH("T")); //0.2);
        public static AAsInfo Valine = new AAsInfo("Valine (V)", "v", getSimulationNEH("V")); //0.56);
        public static AAsInfo Tryptophan = new AAsInfo("Tryptophan (W)", "w", getSimulationNEH("W")); //0.08);
        public static AAsInfo Tyrosine = new AAsInfo("Tyrosine (Y)", "y", getSimulationNEH("y")); //0.42);


        public static List<AAsInfo> getAllAAs()
        {

            List<AAsInfo> list = new List<AAsInfo> { Alanine, Cysteine, Aspartic_acid, Glutamic_acid, Phenylalanine, Glycine, Histidine, Isoleucine, lysine,
                                                        Leucine, Methionine, Asparagine, Proline, Glutamine, Arginine, Serine, Threonine, Valine, Tryptophan, Tyrosine };


            return list;
        }
        #endregion


        public static float[] generateTheroteicalIsotopomerValues(float bwe, float[] fNatIsotopes, float neh)
        {

            float[] fLabIsotopes = new float[6];
            MassIsotopomers MIDyn = new MassIsotopomers();
            MIDyn.CalculateMIDynamics(fNatIsotopes, fLabIsotopes, bwe, neh);
            fLabIsotopes = normalize_list(fLabIsotopes.Select(x => (double)(x)).ToList()).Select(x => (float)x).ToArray();
            return fLabIsotopes;

        }
        public static double computeRsquared(List<double> experimentalValue, List<double> fitvalue)
        {
            var mean_exp = experimentalValue.Where(x => !double.IsNaN(x)).Average();
            double rss = 0;
            double ss = 0;
            double rsquared = double.NaN;

            for (int i = 0; i < experimentalValue.Count; i++)
            {
                if (!double.IsNaN(experimentalValue[i]))
                {
                    ss = ss + Math.Pow((double)(experimentalValue[i] - mean_exp), 2);
                    rss = rss + Math.Pow((double)(experimentalValue[i] - fitvalue[i]), 2);
                }
            }

            //if (r.Rateconst > 0.0006) RSquare = 1 - (rss / ss);
            //else RSquare = 1 - (diff);

            if (ss != 0)
                rsquared = 1 - (rss / ss);

            return rsquared;
        }


        public void compute_two_parm_rate_and_gof_neh(List<DataRecord> data)
        {

            double ph = 1.5574E-4;
            float pw = (float)0.046;
            foreach (var record in data)
            {

                List<double> res = computeTwoCompRate(record.exp_times.ToArray(), record.RIAs.ToArray(), record.M0);
                record.two_param_rate = (float)res[0];
                record.two_param_i0_asymp = (float)res[2];

                List<double> i0_computed = new List<double>();
                for (int i = 0; i < record.exp_times.Count(); i++)
                {
                    var temp_i0_computed = record.two_param_i0_asymp + ((double)record.M0 - record.two_param_i0_asymp) * Math.Exp(-record.two_param_rate * record.exp_times[i]);
                    i0_computed.Add((double)temp_i0_computed);
                }
                record.two_param_rsquared = (float)computeRsquared(record.RIAs.Select(x => (double)x).ToList(), i0_computed);
                //record.two_param_rmse = (float)RMSE(record.RIAs.Select(x => (double)x).ToList(), i0_computed);
                record.NEH_two_param = ((Math.Log((double)record.two_param_i0_asymp / record.M0) / Math.Log(1 - (pw / (1 - ph)))));//(int)Math.Round

                #region compute neh error

                MassIsotopomers MIDyn = new MassIsotopomers();


                float[] fNatIsotopes = new float[6];
                float[] fLabIsotopes = new float[6];
                float[] current_Exp_isotope_distribution = new float[6];

                #region theoretical spec
                fNatIsotopes = new float[6];
                fNatIsotopes = new float[6];
                fNatIsotopes[0] = (float)record.M0 / 100;
                fNatIsotopes[1] = (float)record.M1 / 100;
                fNatIsotopes[2] = (float)record.M2 / 100;
                fNatIsotopes[3] = (float)record.M3 / 100;
                fNatIsotopes[4] = (float)record.M4 / 100;
                fNatIsotopes[5] = (float)record.M5 / 100;
                MIDyn.CalculateMIDynamics(fNatIsotopes, current_Exp_isotope_distribution, record.pxt, (float)record.NEH_sim);
                current_Exp_isotope_distribution = normalize_list(current_Exp_isotope_distribution.Select(x => (double)(x)).ToList()).Select(x => (float)x).ToArray();
                #endregion

                #region experimetnal spec
                fNatIsotopes = new float[6];
                fNatIsotopes = new float[6];
                fNatIsotopes[0] = (float)record.M0 / 100;
                fNatIsotopes[1] = (float)record.M1 / 100;
                fNatIsotopes[2] = (float)record.M2 / 100;
                fNatIsotopes[3] = (float)record.M3 / 100;
                fNatIsotopes[4] = (float)record.M4 / 100;
                fNatIsotopes[5] = (float)record.M5 / 100;
                MIDyn.CalculateMIDynamics(fNatIsotopes, fLabIsotopes, pw, (float)record.NEH_two_param);
                fLabIsotopes = normalize_list(fLabIsotopes.Select(x => (double)(x)).ToList()).Select(x => (float)x).ToArray();
                #endregion

                double[] theoretical_vals_at_t = { fLabIsotopes[0], fLabIsotopes[1], fLabIsotopes[2], fLabIsotopes[3], fLabIsotopes[4], fLabIsotopes[5] };
                var RMSE_error = RMSE(current_Exp_isotope_distribution.Select(x => (double)x).ToList(), theoretical_vals_at_t.ToList());
                record.two_param_rmse = (float)RMSE_error;
                #endregion


            }
        }

        public unsafe static List<double> computeTwoCompRate(float[] TimeCourseDates, float[] TimeCourseI0Isotope, float M0)
        {

            double fDegradationConstant = 0;
            double alpha = 0;
            double io_asy_computed = 0;
            float scale = 60;

            float rkd, rks, fx1;
            float[] scaled_time = TimeCourseDates.Select(x => x / scale).ToArray();

            fixed (float* ptr_TimeCourseDates = scaled_time)
            fixed (float* ptr_TimeCourseI0Isotope = TimeCourseI0Isotope)
            {

                try
                {
                    LBFGS lbfgs = new LBFGS(ptr_TimeCourseDates, TimeCourseDates.Count(), 2, "Two_Parameter_Exponential");
                    lbfgs.InitializeTime();
                    var nRet = lbfgs.Optimize(ptr_TimeCourseI0Isotope, M0, 0, &rkd, &rks, &fx1);

                    fDegradationConstant = Math.Exp(lbfgs.fParams[0]) / scale;
                    alpha = 1 / (1 + Math.Exp(lbfgs.fParams[1]));

                    io_asy_computed = alpha * M0;


                    //Console.WriteLine(TimeCourseI0Isotope[8] + " \t" + io_asy_computed + " \t" + Math.Abs(io_asy_computed - TimeCourseI0Isotope[8]));

                    lbfgs.Release_Memory();
                }
                catch (Exception ex)
                {
                    return null;
                    //continue;
                }


            }

            return new List<double> { fDegradationConstant, alpha, io_asy_computed };
        }



    }


    public class PeptidesPassedNEHfilters
    {

        public string Protein { get; set; }
        public string Peptide { get; set; }
        public int Charge { get; set; }
        public int NumberOfHydrogens { get; set; }
        public double T_NEH { get; set; }
        public double APE_NEH { get; set; }
        public double APE_RMSE { get; set; }
        public double MPE_NEH { get; set; }
        public double MPE_RMSE { get; set; }
        public double Asym_NEH { get; set; }
        public double Asym_RMSE { get; set; }


        public double sd_NEH { get; set; }
        public double sd_RMSE { get; set; }


        public double il_NEH { get; set; }
        public double il_RMSE { get; set; }

        public double i0_0 { get; set; }
        public double i1_0 { get; set; }
        public double i2_0 { get; set; }
        public double i0_31 { get; set; }
        public double i1_31 { get; set; }
        public double i2_31 { get; set; }
    }


    public class AAsInfo
    {
        public AAsInfo(string name, string symbol, double value)
        {
            this.AA_name = name;
            this.AA_symbol = symbol.ToUpper();
            this.AA_NEH_tritium = value;
        }
        public string AA_name { get; set; }
        public string AA_symbol { get; set; }
        public double AA_NEH_tritium { get; set; }
        public double AA_NEH_Asym { get; set; }
        public double AA_NEH_APE { get; set; }
        public double AA_NEH_MPE { get; set; }
        public double AA_NEH_Asym_ci { get; set; }
        public double AA_NEH_APE_ci { get; set; }
        public double AA_NEH_MPE_ci { get; set; }
        public double AA_NEH_APE_MPE { get; set; }
        public double AA_NEH_APE_MPE_ci { get; set; }

    }

    public class DataRecord
    {
        public string Protein { get; set; }
        public string Peptide { get; set; }
        public int Charge { get; set; }
        public double NEH { get; set; }
        public double NEH_sim { get; set; }
        public double NEH_APE { get; set; }
        public double NEH_MPE { get; set; }
        public double NEH_two_param { get; set; }
        public float M0 { get; set; }
        public float M1 { get; set; }
        public float M2 { get; set; }
        public float M3 { get; set; }
        public float M4 { get; set; }
        public float M5 { get; set; }

        public float I0 { get; set; }
        public float I1 { get; set; }
        public float I2 { get; set; }
        public float I3 { get; set; }
        public float I4 { get; set; }
        public float I5 { get; set; }

        public float I0_exp_t { get; set; }
        public float I1_exp_t { get; set; }
        public float I2_exp_t { get; set; }
        public float I3_exp_t { get; set; }
        public float I4_exp_t { get; set; }
        public float I5_exp_t { get; set; }
        public float pxt { get; set; }
        public float two_param_rate { get; set; }
        public float two_param_rsquared { get; set; }
        public float two_param_i0_asymp { get; set; }
        public float two_param_rmse { get; set; }
        public float MPE_rmse { get; set; }
        public float APE_rmse { get; set; }

        public List<float> RIAs = new List<float>();
        public List<float> exp_times = new List<float>();
        public List<List<float>> all_exp_spec = new List<List<float>>();

        public void setExperimentalValues(float[] experimental_data)
        {
            I0_exp_t = experimental_data[0];
            I1_exp_t = experimental_data[1];
            I2_exp_t = experimental_data[2];
            I3_exp_t = experimental_data[3];
            I4_exp_t = experimental_data[4];
            I5_exp_t = experimental_data[5];
        }
    }

}
