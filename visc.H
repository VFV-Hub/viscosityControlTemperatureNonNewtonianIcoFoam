volScalarField viscOld = fluid.nu();
volScalarField viscNew = fluid.nu();

major=0;
viscOld = fluid.nu();

fluid.correct();

viscNew = fluid.nu();
forAll(mesh.cells(),celli)
{
	if(major<=(mag(viscNew[celli]/max(viscOld[celli],VSMALL)-1.)))
	{
		major = mag(viscNew[celli]/max(viscOld[celli],VSMALL)-1.);
	}	
}
