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
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public MRAComponent()
          : base("Multi-body Rope Approach", "MRA",
              "Execute the MRA",
              "MRA", "run_MRA")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Max Lengths", "Max Lengths", "Max Lengths", GH_ParamAccess.list);
            pManager.AddNumberParameter("Ratio Limit", "Ratio Limit", "Ratio Limit", GH_ParamAccess.item);
            pManager.AddNumberParameter("Stiffness", "Stiffness", "Stiffness", GH_ParamAccess.item);
            pManager.AddNumberParameter("Mass", "Mass", "Mass", GH_ParamAccess.item);
            pManager.AddNumberParameter("Time Step", "Time Step", "Time Step", GH_ParamAccess.item);
            pManager.AddNumberParameter("Max Iterations", "Max Iterations", "Max Iterations", GH_ParamAccess.item);
            pManager.AddPointParameter("ConPointS", "ConPointS", "ConPointS", GH_ParamAccess.list);
            pManager.AddPointParameter("ConPointE", "ConPointE", "ConPointE", GH_ParamAccess.list);
            pManager.AddPointParameter("UnConPointS", "UnConPointS", "UnConPointS", GH_ParamAccess.list);
            pManager.AddPointParameter("UnConPointE", "UnConPointE", "UnConPointE", GH_ParamAccess.list);
            pManager.AddPointParameter("Loaded points", "Loaded points", "Loaded points", GH_ParamAccess.list);
            pManager.AddVectorParameter("Load vectors", "Loaded vectors", "Loaded vectors", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stiffness RM", "Stiffness RM", "Stiffness RM", GH_ParamAccess.item);
            pManager.AddNumberParameter("Time Step RM", "Time Step RM", "Time Step RM", GH_ParamAccess.item);
            pManager.AddNumberParameter("Max Iterations RM", "Max Iterations RM", "Max Iterations RM", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Run MRA RM", "Run MRA RM", "Run MRA RM", GH_ParamAccess.item);
            pManager[10].Optional = true;
            pManager[11].Optional = true;
            pManager[12].Optional = true;
            pManager[13].Optional = true;
            pManager[14].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddIntegerParameter("Lasche", "Lasche", "Lasche", GH_ParamAccess.tree);
            pManager.AddIntegerParameter("Num Loose", "Num Loose", "Num Loose", GH_ParamAccess.list);
            pManager.AddNumberParameter("Mean Loose", "Mean Loose", "Mean Loose", GH_ParamAccess.list);
            pManager.AddNumberParameter("Min Loose", "Min Loose", "Min Loose", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Num Over", "Num Over", "Num Over", GH_ParamAccess.list);
            pManager.AddNumberParameter("Mean Over", "Mean Over", "Mean Over", GH_ParamAccess.list);
            pManager.AddNumberParameter("Max Over", "Max Over", "Max Over", GH_ParamAccess.list);
            pManager.AddNumberParameter("Per Dif", "Per Dif", "Per Dif", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Connectivity", "Connectivity", "Connectivity", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Tolerance", "Tolerance", "Tolerance", GH_ParamAccess.item);
            pManager.AddPointParameter("Final Points", "Final Points", "Final Points", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Lfinal", "Lfinal", "Lfinal", GH_ParamAccess.tree);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<double> maxLengths = new List<double>();
            double ratioLimit = new double();
            double stiffness = new double();
            double mass = new double();
            double delta_t = new double();
            double maxIter = new double();
            double stiffnessRM = new double();
            double delta_tRM = new double();
            double maxIterRM = new double();
            bool runMraRM = false;

            DA.GetDataList("Max Lengths", maxLengths);
            DA.GetData("Ratio Limit", ref ratioLimit);
            DA.GetData("Stiffness", ref stiffness);
            DA.GetData("Mass", ref mass);
            DA.GetData("Time Step", ref delta_t);
            DA.GetData("Max Iterations", ref maxIter);
            List<Point3d> iConPointS = new List<Point3d>();
            DA.GetDataList("ConPointS", iConPointS);
            List<Point3d> iConPointE = new List<Point3d>();
            DA.GetDataList("ConPointE", iConPointE);
            List<Point3d> iUnConPointS = new List<Point3d>();
            DA.GetDataList("UnConPointS", iUnConPointS);
            List<Point3d> iUnConPointE = new List<Point3d>();
            DA.GetDataList("UnConPointE", iUnConPointE);
            List<Point3d> loadedPoints = new List<Point3d>();
            DA.GetDataList("Loaded points", loadedPoints);
            List<Vector3d> loadMagnitude = new List<Vector3d>();
            DA.GetDataList("Load vectors", loadMagnitude);
            DA.GetData("Stiffness RM", ref stiffnessRM);
            DA.GetData("Time Step RM", ref delta_tRM);
            DA.GetData("Max Iterations RM", ref maxIterRM);
            DA.GetData("Run MRA RM", ref runMraRM);

            if (runMraRM == true)
            {
                if (stiffnessRM == 0 || delta_tRM == 0 || maxIterRM == 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Assign parameters for the RM method");
                    return;
                }
            }
            if (loadedPoints.Count() != loadMagnitude.Count())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The size of the Loaded Points and Load Magnitude is not the same");
                return;
            }

            // Damping calculation
            damping = 2 * Math.Sqrt(stiffness / mass) * zita;

            // Read Geometry
            List<double[]> start_ext = new List<double[]>();
            List<double[]> end_ext = new List<double[]>();
            List<double[]> start_int = new List<double[]>();
            List<double[]> end_int = new List<double[]>();

            for (int i = 0; i < iConPointS.Count; i++)
            {
                double[] tempArrow1 = new double[3] { iConPointS[i].X, iConPointS[i].Y, iConPointS[i].Z };
                double[] tempArrow2 = new double[3] { iConPointE[i].X, iConPointE[i].Y, iConPointE[i].Z };
                start_ext.Add(tempArrow1);
                end_ext.Add(tempArrow2);
            }
            for (int i = 0; i < iUnConPointS.Count; i++)
            {
                double[] tempArrow3 = new double[3] { iUnConPointS[i].X, iUnConPointS[i].Y, iUnConPointS[i].Z };
                double[] tempArrow4 = new double[3] { iUnConPointE[i].X, iUnConPointE[i].Y, iUnConPointE[i].Z };
                start_int.Add(tempArrow3);
                end_int.Add(tempArrow4);
            }

            // Create external and internal points
            List<double[]> points_ext = Enumerable.Concat(start_ext, end_ext).ToList();
            List<double[]> points_int = Enumerable.Concat(start_int, end_int).ToList();

            // Eliminate doubles in external
            List<double[]> constraints = new List<double[]>();
            var distinct = points_ext.Distinct(new DistinctDoubleArrayComparer());
            foreach (var item in distinct)
            {
                double[] tempArrow = new double[3] { item[0], item[1], item[2] };
                constraints.Add(tempArrow);
            }

            // Eliminate doubles in internal
            List<double[]> unconstraint = new List<double[]>();
            distinct = points_int.Distinct(new DistinctDoubleArrayComparer());
            foreach (var item in distinct)
            {
                double[] tempArrow = new double[3] { item[0], item[1], item[2] };
                unconstraint.Add(tempArrow);
            }

            // Eliminate externals from internals
            List<double[]> tempDuplicates = new List<double[]>();
            var duplicates = constraints.Intersect(unconstraint, new DistinctDoubleArrayComparer());
            foreach (var item in duplicates)
            {
                double[] tempArrow = new double[3] { item[0], item[1], item[2] };
                tempDuplicates.Add(tempArrow);
            }
            foreach (var item in tempDuplicates)
            {
                unconstraint.RemoveAll(arr => arr.SequenceEqual(item));
            }
            List<double[]> points = Enumerable.Concat(constraints, unconstraint).ToList();
            n_point = points.Count();
            n_vinc = constraints.Count();
            n_no_vinc = unconstraint.Count();
            n_ropes = start_int.Count();

            // Connectivity - typical and advanced formation
            List<int[]> link = new List<int[]>();
            for (int i = 0; i < n_ropes; i++)
            {
                double[] tempArrow = start_int[i];
                int index1 = points.FindIndex(l => Enumerable.SequenceEqual(tempArrow, l));
                tempArrow = end_int[i];
                int index2 = points.FindIndex(l => Enumerable.SequenceEqual(tempArrow, l));
                int[] tempIntArrow = new int[2] { index1, index2 };
                link.Add(tempIntArrow);
            }
            List<int[]> link_p = new List<int[]>();
            AdvancedConnectivity(link, ref link_p, n_point, n_ropes);

            // Gravity Loads
            List<double[]> carichi = new List<double[]>();
            for (int i = 0; i < n_point; i++)
            {
                double[] tempCarichi;
                if (i < n_vinc)
                {
                    tempCarichi = new double[3] { 0, 0, 0 };
                    carichi.Add(tempCarichi);
                }
                else
                {
                    tempCarichi = new double[3] { 0, 0, mass * g };
                    carichi.Add(tempCarichi);
                }
            }

            // Point Loads
            int tempId;
            for (int i = 0; i < loadedPoints.Count; i++)
            {
                double[] tempArrow1 = new double[3] { loadedPoints[i].X, loadedPoints[i].Y, loadedPoints[i].Z };
                tempId = points.FindIndex(l => Enumerable.SequenceEqual(tempArrow1, l));
                double[] tempArrow2 = new double[3] { loadMagnitude[i].X, loadMagnitude[i].Y, loadMagnitude[i].Z };
                for (int j = 0; j < 3; j++)
                {
                    carichi[tempId][j] += tempArrow2[j];
                }
            }

            // Tranform from List to arrays to continue with the iterative method
            double[,] pointsArray = new double[n_point, 3];
            for (int i = 0; i < n_point; i++)
            {
                pointsArray[i, 0] = points[i][0];
                pointsArray[i, 1] = points[i][1];
                pointsArray[i, 2] = points[i][2];
            }

            // Initialization velocity, acceleration, lengths
            DataTempLoop myDataTemp = new DataTempLoop();
            myDataTemp.lCurrent = RopeLength(link, pointsArray, n_ropes);
            myDataTemp.sCurrent = new double[n_ropes];
            myDataTemp.uCurrent = MyDeepCopyArrayDouble(pointsArray, n_point);
            myDataTemp.vCurrent = new double[n_point, 3];
            myDataTemp.aCurrent = new double[n_point, 3];
            myDataTemp.lAfter = new double[n_ropes];
            myDataTemp.sAfter = new double[n_ropes];
            myDataTemp.uAfter = new double[n_point, 3];
            myDataTemp.vAfter = new double[n_point, 3];
            myDataTemp.aAfter = new double[n_point, 3];
            myDataTemp.norm_v = new double[n_point];
            myDataTemp.norm_a = new double[n_point];

            // N order of rope - limit definition
            total_orders = maxLengths.Count();
            List<double> l_limit_total = new List<double>();
            if (total_orders > 1)
            {
                for (int i = 1; i < total_orders; i++)
                {
                    l_limit_total.Add(Math.Round(ratioLimit * (maxLengths[i - 1] - maxLengths[i]) + maxLengths[i], 2));
                }
            }

            // Start - Iterating loop
            List<FinalResults> myFinalResults = new List<FinalResults>();
            List<DataTempLoop> myAllDataTemp = new List<DataTempLoop>();
            for (int temp_order = 0; temp_order < total_orders; temp_order++)
            {
                List<double> l_limit_n = new List<double>();
                if (temp_order > 0)
                {
                    myDataTemp = new DataTempLoop();
                    myDataTemp.lCurrent = MyDeepCopyDouble(myAllDataTemp[temp_order - 1].lCurrent);
                    myDataTemp.sCurrent = new double[n_ropes];
                    myDataTemp.uCurrent = MyDeepCopyArrayDouble(myAllDataTemp[temp_order - 1].uCurrent, n_point);
                    myDataTemp.vCurrent = new double[n_point, 3];
                    myDataTemp.aCurrent = new double[n_point, 3];
                    myDataTemp.lAfter = new double[n_ropes];
                    myDataTemp.sAfter = new double[n_ropes];
                    myDataTemp.uAfter = new double[n_point, 3];
                    myDataTemp.vAfter = new double[n_point, 3];
                    myDataTemp.aAfter = new double[n_point, 3];
                    myDataTemp.norm_v = new double[n_point];
                    myDataTemp.norm_a = new double[n_point];
                    l_limit_n = l_limit_total.GetRange(0, temp_order);
                }

                int iter = 1;
                // Main loop
                while (iter < maxIter)
                {
                    iter += 1;

                    // Calculate forces
                    double[] s_temp = CalcForce(n_ropes, temp_order, myDataTemp.lCurrent, stiffness, maxLengths, l_limit_n);
                    myDataTemp.sAfter = MyDeepCopyDouble(s_temp);

                    // Free masses
                    FreeMass(ref myDataTemp, n_vinc, n_point, carichi, mass, link_p, link, myDataTemp.lCurrent, damping, delta_t);

                    // Constraint masses
                    for (int i = 0; i < n_vinc; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            double tempValue = myDataTemp.uCurrent[i, j];
                            myDataTemp.uAfter[i, j] = tempValue;
                        }
                    }

                    // New Length
                    myDataTemp.lAfter = RopeLength(link, myDataTemp.uAfter, n_ropes);

                    // Tollerance check
                    double tempNormV = NormEvaluation(myDataTemp.norm_v);
                    double tempNormA = NormEvaluation(myDataTemp.norm_a);
                    if (tempNormV < tollLoop && tempNormA < tollLoop)
                    {
                        break;
                    }

                    // End of while loop
                    myDataTemp.lCurrent = MyDeepCopyDouble(myDataTemp.lAfter);
                    myDataTemp.sCurrent = MyDeepCopyDouble(myDataTemp.sAfter);
                    myDataTemp.uCurrent = MyDeepCopyArrayDouble(myDataTemp.uAfter, n_point);
                    myDataTemp.vCurrent = MyDeepCopyArrayDouble(myDataTemp.vAfter, n_point);
                    myDataTemp.aCurrent = MyDeepCopyArrayDouble(myDataTemp.aAfter, n_point);
                }
                FinalResults tempFinalResults = LengthControl(n_ropes, temp_order, myDataTemp.lAfter, maxLengths, toll_c);
                tempFinalResults.uFinal = MyDeepCopyArrayDouble(myDataTemp.uCurrent, n_point);
                myFinalResults.Add(tempFinalResults);
                myAllDataTemp.Add(myDataTemp);
            }

            int total_orders_rm = total_orders;
            if (runMraRM == true)
            {
                List<double> l_limit_n = new List<double>();
                myDataTemp = new DataTempLoop();
                myDataTemp.lCurrent = MyDeepCopyDouble(myAllDataTemp[total_orders - 1].lCurrent);
                myDataTemp.sCurrent = new double[n_ropes];
                myDataTemp.uCurrent = MyDeepCopyArrayDouble(myAllDataTemp[total_orders - 1].uCurrent, n_point);
                myDataTemp.vCurrent = new double[n_point, 3];
                myDataTemp.aCurrent = new double[n_point, 3];
                myDataTemp.lAfter = new double[n_ropes];
                myDataTemp.sAfter = new double[n_ropes];
                myDataTemp.uAfter = new double[n_point, 3];
                myDataTemp.vAfter = new double[n_point, 3];
                myDataTemp.aAfter = new double[n_point, 3];
                myDataTemp.norm_v = new double[n_point];
                myDataTemp.norm_a = new double[n_point];
                l_limit_n = l_limit_total.GetRange(0, total_orders - 1);

                int iter = 1;
                // Main loop
                while (iter < maxIterRM)
                {
                    iter += 1;

                    // Calculate forces
                    double[] s_temp = CalcForceRM(n_ropes, total_orders - 1, myDataTemp.lCurrent, stiffness, stiffnessRM, maxLengths, l_limit_n, toll_c);
                    myDataTemp.sAfter = MyDeepCopyDouble(s_temp);

                    // Free masses
                    FreeMass(ref myDataTemp, n_vinc, n_point, carichi, mass, link_p, link, myDataTemp.lCurrent, damping, delta_tRM);

                    // Constraint masses
                    for (int i = 0; i < n_vinc; i++)
                    {
                        for (int j = 0; j < 3; j++)
                        {
                            double tempValue = myDataTemp.uCurrent[i, j];
                            myDataTemp.uAfter[i, j] = tempValue;
                        }
                    }

                    // New Length
                    myDataTemp.lAfter = RopeLength(link, myDataTemp.uAfter, n_ropes);

                    // Tollerance check
                    double tempNormV = NormEvaluation(myDataTemp.norm_v);
                    double tempNormA = NormEvaluation(myDataTemp.norm_a);
                    if (tempNormV < tollLoop && tempNormA < tollLoop)
                    {
                        break;
                    }

                    // End of while loop
                    myDataTemp.lCurrent = MyDeepCopyDouble(myDataTemp.lAfter);
                    myDataTemp.sCurrent = MyDeepCopyDouble(myDataTemp.sAfter);
                    myDataTemp.uCurrent = MyDeepCopyArrayDouble(myDataTemp.uAfter, n_point);
                    myDataTemp.vCurrent = MyDeepCopyArrayDouble(myDataTemp.vAfter, n_point);
                    myDataTemp.aCurrent = MyDeepCopyArrayDouble(myDataTemp.aAfter, n_point);
                }
                FinalResults tempFinalResults = LengthControl(n_ropes, total_orders - 1, myDataTemp.lAfter, maxLengths, toll_c);
                tempFinalResults.uFinal = MyDeepCopyArrayDouble(myDataTemp.uCurrent, n_point);
                myFinalResults.Add(tempFinalResults);
                myAllDataTemp.Add(myDataTemp);
                total_orders_rm += 1;
            }

            //Assign to Grasshoper - Output
            List<int> numLooseList = new List<int>();
            List<double> meanLooseList = new List<double>();
            List<double> minLooseList = new List<double>();
            List<int> numOverList = new List<int>();
            List<double> meanOverList = new List<double>();
            List<double> maxOverList = new List<double>();
            List<double> perDiffList = new List<double>();
            for (int i = 0; i < total_orders_rm; i++)
            {
                numLooseList.Add(myFinalResults[i].numLoose);
                meanLooseList.Add(myFinalResults[i].meanLoose);
                minLooseList.Add(myFinalResults[i].minLoose);
                numOverList.Add(myFinalResults[i].numOver);
                meanOverList.Add(myFinalResults[i].meanOver);
                maxOverList.Add(myFinalResults[i].maxOver);
                perDiffList.Add(myFinalResults[i].perDiff);
            }

            DA.SetDataList(1, numLooseList);
            DA.SetDataList(2, meanLooseList);
            DA.SetDataList(3, minLooseList);
            DA.SetDataList(4, numOverList);
            DA.SetDataList(5, meanOverList);
            DA.SetDataList(6, maxOverList);
            DA.SetDataList(7, perDiffList);

            Grasshopper.DataTree<int> myLasche = new Grasshopper.DataTree<int>();
            for (int i = 0; i < total_orders_rm; i++)
            {
                List<int> temp = myFinalResults[i].lasche1;
                Grasshopper.Kernel.Data.GH_Path pathi = new Grasshopper.Kernel.Data.GH_Path(i);
                myLasche.AddRange(temp, pathi);
            }
            DA.SetDataTree(0, myLasche);

            // Collect info for plotting
            Grasshopper.DataTree<int> myConnectivity = new Grasshopper.DataTree<int>();
            for (int i = 0; i < n_ropes; i++)
            {
                List<int> temp = new List<int>();
                temp.Add(link[i][0]);
                temp.Add(link[i][1]);
                Grasshopper.Kernel.Data.GH_Path pathi = new Grasshopper.Kernel.Data.GH_Path(i);
                myConnectivity.AddRange(temp, pathi);
            }
            DA.SetDataTree(8, myConnectivity);
            DA.SetData(9, toll_c);

            Grasshopper.DataTree<Point3d> finalPoints = new Grasshopper.DataTree<Point3d>();
            for (int j = 0; j < myFinalResults.Count; j++)
            {
                List<Point3d> finalPoints1 = new List<Point3d>();
                double[,] uFinal1 = myFinalResults[j].uFinal;
                for (int i = 0; i < n_point; i++)
                {
                    Point3d tempPoint = new Point3d();
                    tempPoint.X = uFinal1[i, 0];
                    tempPoint.Y = uFinal1[i, 1];
                    tempPoint.Z = uFinal1[i, 2];
                    finalPoints1.Add(tempPoint);
                }
                Grasshopper.Kernel.Data.GH_Path pathi = new Grasshopper.Kernel.Data.GH_Path(j);
                finalPoints.AddRange(finalPoints1, pathi);
            }
            DA.SetDataTree(10, finalPoints);

            Grasshopper.DataTree<double> finalLengths = new Grasshopper.DataTree<double>();
            for (int j = 0; j < myFinalResults.Count; j++)
            {
                List<double> finalLengths1 = myFinalResults[j].lFinal;
                Grasshopper.Kernel.Data.GH_Path pathi = new Grasshopper.Kernel.Data.GH_Path(j);
                finalLengths.AddRange(finalLengths1, pathi);
            }
            DA.SetDataTree(11, finalLengths);

        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                // You can add image files to your project resources and access them like this:
                //return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6122ee31-5117-4f92-9632-3775d47ddbb0"); }
        }
    }
}
