using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;


// In order to load the result of this wizard, you will also need to
// add the output bin/ folder of this project to the list of loaded
// folder in Grasshopper.
// You can use the _GrasshopperDeveloperSettings Rhino command for that.

namespace MRA
{
    public partial class MRAComponent : GH_Component
    {
        protected double g = 9.81;        // Gravity acceleration g - m/s²
        protected double zita = 0.95;     // 1/s - Coefficient damping/mass
        protected double toll_c = 0.01;   // 2-5% of the rope's length
        protected double tollLoop = 0.001; // Tollerance for finishing the while loop
        protected double damping;         // Damping
        protected int n_point;            // Number of points
        protected int n_vinc;             // Number of constrained points
        protected int n_no_vinc;          // Number of non constrained points
        protected int n_ropes;            // Number of total ropes
        protected int total_orders;       // Number of total orders

        /// <summary>
        /// Create advanced connectivity
        /// </summary>
        /// <param name="link1"><Typical connectivity>
        /// <param name="link_p1"><Advanced connectivity>
        protected void AdvancedConnectivity(List<int[]> link1, ref List<int[]> link_p1, int n_point1, int n_ropes1)
        {
            for (int i = 0; i < n_point1; i++)
            {
                int[] tempArrow1 = new int[6] { -1, -1, -1, -1, -1, -1 };
                int k = 0;
                for (int j = 0; j < n_ropes1; j++)
                {
                    int[] tempArrow2 = link1[j];
                    if (i == tempArrow2[0])
                    {
                        tempArrow1[k] = tempArrow2[1];
                        k += 1;
                    }
                    else if (i == tempArrow2[1])
                    {
                        tempArrow1[k] = tempArrow2[0];
                        k += 1;
                    }
                }
                link_p1.Add(tempArrow1);
            }
        }

        /// <summary>
        /// Evaluate ropes' lengths
        /// </summary>
        /// <param name="link1"></param>
        /// <param name="points1"></param>
        /// <param name="n_ropes1"></param>
        /// <returns></returns>
        protected double[] RopeLength(List<int[]> link1, double[,] points1, int n_ropes1)
        {
            double[] lengths1 = new double[n_ropes1];
            double tempLength;
            for (int i = 0; i < n_ropes1; i++)
            {
                int[] temp1 = link1[i];
                double x1 = points1[temp1[0], 0];
                double y1 = points1[temp1[0], 1];
                double z1 = points1[temp1[0], 2];
                double x2 = points1[temp1[1], 0];
                double y2 = points1[temp1[1], 1];
                double z2 = points1[temp1[1], 2];
                tempLength = Math.Sqrt(Math.Pow((x2 - x1), 2) + Math.Pow((y2 - y1), 2) + Math.Pow((z2 - z1), 2));
                lengths1[i] = tempLength;
            }

            return lengths1;
        }

        /// <summary>
        /// Calculate forces in categories for n - orders
        /// </summary>
        /// <param name="n_ropes1"></param>
        /// <param name="n_order"></param>
        /// <param name="l_temp1"></param>
        /// <param name="stiffness1"></param>
        /// <param name="maxLengths1"></param>
        /// <param name="l_limit_total1"></param>
        /// <returns></returns>
        protected double[] CalcForce(int n_ropes1, int n_order, double[] l_temp1, double stiffness1, List<double> maxLengths1, List<double> l_limit_total1)
        {
            double[] s = new double[n_ropes1];
            for (int j = 0; j < n_ropes1; j++)
            {
                s[j] = 0.0;
                if (l_temp1[j] > maxLengths1[0])
                {
                    s[j] = Math.Round(stiffness1 * (l_temp1[j] - maxLengths1[0]), 4);
                    continue;
                }
                if (n_order > 0)
                {
                    double lmax_n;
                    double l_limit_n;
                    for (int i = 1; i <= n_order; i++)
                    {
                        lmax_n = maxLengths1[i];
                        l_limit_n = l_limit_total1[i - 1];
                        if (l_temp1[j] > lmax_n && l_temp1[j] < l_limit_n)
                        {
                            s[j] = Math.Round(stiffness1 * (l_temp1[j] - lmax_n), 4);
                            continue;
                        }
                    }
                }
            }

            return s;
        }

        /// <summary>
        /// Calculate forces for repulsive masses
        /// </summary>
        protected double[] CalcForceRM(int n_ropes1, int n_order, double[] l_temp1, double stiffness1, double stiffness2, List<double> maxLengths1, List<double> l_limit_total1, double toll_c1)
        {
            double[] s = new double[n_ropes1];
            for (int j = 0; j < n_ropes1; j++)
            {
                s[j] = 0.0;
                if (l_temp1[j] > maxLengths1[0])
                {
                    s[j] = Math.Round(stiffness1 * (l_temp1[j] - maxLengths1[0]), 4);
                    continue;
                }
                double lmax_n;
                double l_limit_n;
                for (int i = 1; i <= n_order; i++)
                {
                    lmax_n = maxLengths1[i];
                    l_limit_n = l_limit_total1[i - 1];
                    if (l_temp1[j] > lmax_n && l_temp1[j] < l_limit_n)
                    {
                        s[j] = Math.Round(stiffness1 * (l_temp1[j] - lmax_n), 4);
                        continue;
                    }
                }

                // Repulsive
                if (l_temp1[j] < maxLengths1[n_order] - toll_c1)
                {
                    s[j] = Math.Round(stiffness2 * (maxLengths1[n_order] - l_temp1[j]), 4);
                    continue;
                }
                for (int i = 1; i <= n_order; i++)
                {
                    lmax_n = maxLengths1[i - 1];
                    l_limit_n = l_limit_total1[i - 1];
                    if (l_temp1[j] > l_limit_n && l_temp1[j] < lmax_n - toll_c1)
                    {
                        s[j] = Math.Round(stiffness2 * (lmax_n - l_temp1[j]), 4);
                        continue;
                    }
                }
            }

            return s;
        }

        /// <summary>
        /// Evaluate new dsiplacements, velocity, accelerations
        /// </summary>
        /// <param name="myDataTemp1"></param>
        /// <param name="n_vinc1"><Number of constraint points>
        /// <param name="n_point1"><Number of points>
        /// <param name="carichi1"><Vector of loads>
        /// <param name="mass1"></Mass>
        /// <param name="link_p1"><Advanced connectivity>
        /// <param name="link1"><Typical connectivity>
        /// <param name="l_temp1"><Vector of current lengths>
        /// <param name="damping1"></Damping>
        /// <param name="delta_t1"><Time step>
        protected void FreeMass(ref DataTempLoop myDataTemp1, int n_vinc1, int n_point1, List<double[]> carichi1, double mass1, List<int[]> link_p1, List<int[]> link1, double[] l_temp1, double damping1, double delta_t1)
        {
            double[] p_i;
            double[] u_i = new double[3];
            double[] v_i = new double[3];
            double[] cost;
            int xA;
            double[] u_A = new double[3];
            double l_A = 0;
            double s_A = 0; ;
            myDataTemp1.uAfter = new double[n_point1, 3];
            myDataTemp1.vAfter = new double[n_point1, 3];
            myDataTemp1.aAfter = new double[n_point1, 3];
            for (int i = n_vinc1; i < n_point1; i++)
            {
                p_i = carichi1[i];
                for (int ijk = 0; ijk < 3; ijk++)
                {
                    u_i[ijk] = myDataTemp1.uCurrent[i, ijk];
                    v_i[ijk] = myDataTemp1.vCurrent[i, ijk];
                }
                cost = MatrixDivisionWithNumber(p_i, mass1);
                // Identify the nodes connected to the i node and concequently calculate the applied force
                for (int j = 0; j < link_p1.Count; j++)
                {
                    xA = link_p1[i][j];
                    if (xA == -1)
                    {
                        break;
                    }
                    for (int ijk = 0; ijk < 3; ijk++)
                    {
                        u_A[ijk] = myDataTemp1.uCurrent[xA, ijk];
                    }
                    for (int k = 0; k < link1.Count; k++)
                    {
                        if (link1[k][0] == i && link1[k][1] == xA)
                        {
                            l_A = l_temp1[k];
                            s_A = myDataTemp1.sAfter[k];
                            break;
                        }
                        else if (link1[k][1] == i && link1[k][0] == xA)
                        {
                            l_A = l_temp1[k];
                            s_A = myDataTemp1.sAfter[k];
                            break;
                        }
                    }
                    for (int ij = 0; ij < 3; ij++)
                    {
                        cost[ij] = cost[ij] - 1 / mass1 * ((u_i[ij] - u_A[ij]) / l_A) * s_A;
                    }
                }

                // Equation: a(i) + zita*v(i) = cost
                // Position at i
                double temp_u1;
                double temp_u2;
                double temp_u3;
                double[] tempNormV = new double[3];
                double[] tempNormA = new double[3];
                for (int j = 0; j < 3; j++)
                {
                    temp_u1 = (Math.Pow(damping1, 2) * u_i[j] + damping1 * v_i[j] - cost[j]) / Math.Pow(damping1, 2);
                    temp_u2 = (Math.Exp(-damping1 * delta_t1) * (damping1 * v_i[j] - cost[j])) / Math.Pow(damping1, 2);
                    temp_u3 = (cost[j] * delta_t1) / damping1;
                    myDataTemp1.uAfter[i, j] = temp_u1 - temp_u2 + temp_u3;
                    myDataTemp1.vAfter[i, j] = (temp_u1 - temp_u2 + temp_u3 - u_i[j]) / delta_t1;
                    myDataTemp1.aAfter[i, j] = (((temp_u1 - temp_u2 + temp_u3 - u_i[j]) / delta_t1) - v_i[j]) / delta_t1;
                    tempNormV[j] = myDataTemp1.vAfter[i, j];
                    tempNormA[j] = myDataTemp1.aAfter[i, j];
                }
                myDataTemp1.norm_v[i] = NormEvaluation(tempNormV);
                myDataTemp1.norm_a[i] = NormEvaluation(tempNormA);
            }
            return;
        }

        /// <summary>
        /// Divide ropes in categories for n - orders
        /// </summary>
        /// <param name="n_ropes1"></param>
        /// <param name="n_order"></param>
        /// <param name="l_temp1"></param>
        /// <param name="maxLengths1"></param>
        /// <param name="toll_c1"></param>
        /// <param name="link1"></param>
        /// <returns></returns>
        protected FinalResults LengthControl(int n_ropes1, int n_order, double[] l_temp1, List<double> maxLengths1, double toll_c1)
        {
            FinalResults finalResults = new FinalResults();
            List<double> loose = new List<double>();
            List<int> lasche = new List<int>();
            List<double> over = new List<double>();
            finalResults.lFinal = new List<double>();

            for (int j = 0; j < n_ropes1; j++)
            {
                if (Math.Round(l_temp1[j], 2) > Math.Round(maxLengths1[0] + toll_c1, 2))
                {
                    over.Add(l_temp1[j]);
                    continue;
                }
                if (Math.Round(l_temp1[j], 2) < Math.Round(maxLengths1[n_order] - toll_c1, 2))
                {
                    loose.Add(l_temp1[j]);
                    lasche.Add(j);
                    continue;
                }
                if (n_order > 0)
                {
                    for (int i = 0; i < n_order; i++)
                    {
                        double lmax_n = maxLengths1[i + 1];
                        double lmax_n_1 = maxLengths1[i];
                        if (Math.Round(l_temp1[j], 2) < Math.Round(lmax_n_1 - toll_c1, 2) && Math.Round(l_temp1[j], 2) > Math.Round(lmax_n + toll_c1, 2))
                        {
                            loose.Add(l_temp1[j]);
                            lasche.Add(j);
                            continue;
                        }
                    }
                }
            }

            finalResults.lasche1 = lasche;
            finalResults.numLoose = loose.Count();
            if (loose.Count() > 0)
            {
                finalResults.meanLoose = loose.Average();
                finalResults.minLoose = loose.Min();
            }
            finalResults.numOver = over.Count();
            if (over.Count() > 0)
            {
                finalResults.meanOver = over.Average();
                finalResults.maxOver = over.Max();
            }
            finalResults.perDiff = 100 * (finalResults.numLoose + finalResults.numOver) / n_ropes1;

            for (int i = 0; i < l_temp1.Length; i++)
            {
                double tempVal = l_temp1[i];
                finalResults.lFinal.Add(tempVal);
            }

            return finalResults;
        }

        /// <summary>
        /// Divide a vector with a number
        /// </summary>
        /// <param name="initMatrix"></param>
        /// <param name="myNumber"></param>
        /// <returns></returns>
        public double[] MatrixDivisionWithNumber(double[] initMatrix, double myNumber)
        {
            double[] finalMatrix = new double[initMatrix.Length];
            for (int i = 0; i < initMatrix.Length; i++)
            {
                finalMatrix[i] = initMatrix[i] / myNumber;
            }

            return finalMatrix;
        }

        /// <summary>
        /// Copy a vector
        /// </summary>
        /// <param name="myList"></param>
        /// <returns></returns>
        protected double[] MyDeepCopyDouble(double[] myList)
        {
            double[] myFinalList = new double[myList.Length];
            for (int i = 0; i < myList.Length; i++)
            {
                double tempValue = myList[i];
                myFinalList[i] = tempValue;
            }

            return myFinalList;
        }

        /// <summary>
        /// Copy a matrix num_row * 3
        /// </summary>
        /// <param name="myList"></param>
        /// <param name="num_row"></param>
        /// <returns></returns>
        protected double[,] MyDeepCopyArrayDouble(double[,] myList, int num_row)
        {
            double[,] myFinalList = new double[num_row, 3];
            for (int i = 0; i < num_row; i++)
            {
                double[] tempValue = new double[3] { myList[i, 0], myList[i, 1], myList[i, 2] };
                myFinalList[i, 0] = tempValue[0];
                myFinalList[i, 1] = tempValue[1];
                myFinalList[i, 2] = tempValue[2];
            }

            return myFinalList;
        }

        /// <summary>
        /// Evaluate the magnitude of a vector
        /// </summary>
        /// <param name="myVector"></param>
        /// <returns></returns>
        protected double NormEvaluation(double[] myVector)
        {
            double norm_v1;
            double temp = 0;
            for (int i = 0; i < myVector.Length; i++)
            {
                temp += Math.Pow(myVector[i], 2);
            }
            norm_v1 = Math.Sqrt(temp);

            return norm_v1;
        }
    }
}
