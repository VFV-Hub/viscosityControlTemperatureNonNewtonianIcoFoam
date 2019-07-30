volScalarField viscOld = fluid.nu();
volScalarField viscNew = fluid.nu();

max=0;
viscOld = fluid.nu();

fluid.correct();

viscNew = fluid.nu();
forAll(mesh.cells(),celli)
{
	if(max<=(mag(viscNew[celli]-viscOld[celli])))
	{
		max = mag(viscNew[celli]-viscOld[celli]);
	}
	//necessário para o primeiro passo de tempo
	if(max==0)
		max=1;		
}
