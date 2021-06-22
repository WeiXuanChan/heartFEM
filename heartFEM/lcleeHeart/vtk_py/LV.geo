//Geometry.HideCompounds = 1;

Mesh.CharacteristicLengthFactor = <<mesh_d>>;  // Coarse

Mesh.Algorithm    = 2; // (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad) (Default=2)
Mesh.RecombineAll = 0;

//Mesh.RemeshAlgorithm = 1; // (0=no split, 1=automatic, 2=automatic only with metis) (Default=0)

//Mesh.RemeshParametrization = 7; // (0=harmonic_circle, 1=conformal_spectral, 2=rbf, 3=harmonic_plane, 4=convex_circle, 5=convex_plane, 6=harmonic square, 7=conformal_fe) (Default=4)

Mesh.Algorithm3D    = 4; // (1=Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree) (Default=1)
Mesh.Recombine3DAll = 0;

Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;

Mesh.Smoothing = 0;

Merge "<<Endofilename>>";
Merge "<<Epifilename>>";

DefineConstant[
  // Angle between two triangles above which an edge is considered as sharp
  angle = {40, Min 20, Max 120, Step 1,
    Name "Parameters/Angle for surface detection"},
  // For complex geometries, patches can be too complex, too elongated or too
  // large to be parametrized; setting the following option will force the
  // creation of patches that are amenable to reparametrization:
  forceParametrizablePatches = {0, Choices{0,1},
    Name "Parameters/Create surfaces guaranteed to be parametrizable"},
  // For open surfaces include the boundary edges in the classification process:
  includeBoundary = 1,
  // Force curves to be split on given angle:
  curveAngle = 180
];
ClassifySurfaces{angle * Pi/180, includeBoundary, forceParametrizablePatches,
                 curveAngle * Pi / 180};
CreateGeometry;

CreateTopology;

ss[] = Surface "*";

L_LV_base = newll; Curve Loop(L_LV_base) = {Boundary{ Surface {ss[0]}; }};

L_epi_base = newll; Curve Loop(L_epi_base) = {Boundary{ Surface {ss[1]}; }};

Physical Surface("LV") = {ss[0]};
Physical Surface("epi") = {ss[1]};

S_base = news; Plane Surface(S_base) = { L_epi_base,L_LV_base};

Physical Surface("BASE") = {S_base};

SL_wall = newsl; Surface Loop(SL_wall) = {ss[1],ss[0], S_base};

V_wall = newv; Volume(V_wall) = {SL_wall};

Physical Volume("WALL") = {V_wall};


//
//// P1 = newp; Point(P1) = {36., 74., 70.};
//// Field[1] = Attractor;
//// Field[1].NodesList = {P1};
//// Field[2] = Threshold;
//// Field[2].IField = 1;
//// Field[2].LcMin = 100.;
//// Field[2].LcMax = 100.;
//// Field[2].DistMin = 0.;
//// Field[2].DistMax = 10.;
//// Background Field = 2;
//
//// Field[1] = Box;
//// Field[1].VIn = 5.;
//// Field[1].VOut = 5.;
//// Field[1].XMin = 30.; 
//// Field[1].XMax = 40.;
//// Field[1].YMin = 70.;
//// Field[1].YMax = 80.;
//// Field[1].ZMin = 69.;
//// Field[1].ZMax = 71.;
//// Background Field = 1;








//+
SetFactory("Built-in");
//+
SetFactory("Built-in");
