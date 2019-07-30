// Momentum predictor
fvVectorMatrix UEqn
(
	fvm::ddt(U)
  + fvm::div(phi, U)
  - fvm::laplacian(fluid.nu(), U)
  - (fvc::grad(U) & fvc::grad(fluid.nu()))
);

if (piso.momentumPredictor())
{
	solve(UEqn == -fvc::grad(p));
}


