using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using System;
using System.Collections.Generic;

namespace MRA
{
    public class GhcMRAplot : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the GhcMRAplot class.
        /// </summary>
        public GhcMRAplot()
          : base("GhcMRAplot", "MRAplot",
              "Plot MRA results",
              "MRA", "run_MRA")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Max Lengths", "Max Lengths", "Max Lengths", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Display Order", "Display Order", "Display Order", GH_ParamAccess.item);
            pManager.AddNumberParameter("Tollerance", "Tollerance", "Tollerance", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Connectivity", "Connectivity", "Connectivity", GH_ParamAccess.tree);
            pManager.AddPointParameter("Final Points", "Final Points", "Final Points", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Lfinal", "Lfinal", "Lfinal", GH_ParamAccess.tree);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Order Points", "Order Points", "Order Points", GH_ParamAccess.list);
            pManager.AddLineParameter("Ropes", "Ropes", "Ropes", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<double> maxLengths = new List<double>();
            int displayOrder = new int();
            double toll_c = new double();
            Grasshopper.Kernel.Data.GH_Structure<GH_Integer> myConnectivity;
            Grasshopper.Kernel.Data.GH_Structure<GH_Point> finalPoints;
            Grasshopper.Kernel.Data.GH_Structure<GH_Number> finalLengths;


            DA.GetDataList(0, maxLengths);
            DA.GetData(1, ref displayOrder);
            DA.GetData(2, ref toll_c);
            DA.GetDataTree(3, out myConnectivity);
            DA.GetDataTree(4, out finalPoints);
            DA.GetDataTree(5, out finalLengths);

            // Retrieve info for plotting
            int n_ropes = myConnectivity.Branches.Count;
            int[,] myConn = new int[n_ropes, 2];
            for (int i = 0; i < n_ropes; i++)
            {
                GH_Integer[] mytemp = myConnectivity.Branches[i].ToArray();
                myConn[i, 0] = mytemp[0].Value;
                myConn[i, 1] = mytemp[1].Value;
            }

            GH_Point[] currentPoints = finalPoints.Branches[displayOrder].ToArray();
            GH_Number[] currentLengths = finalLengths.Branches[displayOrder].ToArray();

            List<Point3d> currentPointsL = new List<Point3d>();
            for (int i = 0; i < currentPoints.Length; i++)
            {
                currentPointsL.Add(currentPoints[i].Value);
            }
            DA.SetDataList(0, currentPointsL);

            List<Line> currentFrames = new List<Line>();
            for (int i = 0; i < n_ropes; i++)
            {
                Line tempLine = new Line(currentPointsL[myConn[i, 0]], currentPointsL[myConn[i, 1]]);
                currentFrames.Add(tempLine);
            }
            DA.SetDataList(1, currentFrames);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("b5468dd5-fd4e-4ea4-bc42-09abdef43a9c"); }
        }
    }
}