from dolfin import *
import logging
logger = logging.getLogger(__name__)


class NSolver(object):

    def __init__(self, params):

        self.parameters = params
        self.isfirstiteration = 0

    def default_parameters(self):
        return {"rel_tol": 1e-6,
                "abs_tol": 1e-6,
                "max_iter": 50}

    def solvenonlinear(self):

        abs_tol = self.default_parameters()["abs_tol"]
        rel_tol = self.default_parameters()["rel_tol"]
        maxiter = self.default_parameters()["max_iter"]
        Jac = self.parameters["Jacobian"]
        Ftotal = self.parameters["F"]
        w = self.parameters["w"]
        bcs = self.parameters["boundary_conditions"]
        solvertype = self.parameters["Type"]

        mesh = self.parameters["mesh"]
        comm = w.function_space().mesh().mpi_comm()

        if(solvertype == 0):
            solve(Ftotal == 0, w, bcs, J=Jac, solver_parameters={"newton_solver": {
                  "relative_tolerance": 1e-6, "absolute_tolerance": 1e-6, "maximum_iterations": maxiter}})
        else:

            it = 0
            if(self.isfirstiteration == 0):
                A, b = assemble_system(Jac, -Ftotal, bcs)
                resid0 = b.norm("l2")
                rel_res = b.norm("l2")/resid0
                res = resid0
                if(MPI.rank(comm) == 0):
                    logger.debug("Iteration: %d, Residual: %.3e, Relative residual: %.3e" % (
                        it, res, rel_res))
                solve(A, w.vector(), b)
                it += 1

            self.isfirstiteration = 1

            B = assemble(Ftotal)
            for bc in bcs:
                bc.apply(B)

            rel_res = 1.0
            res = B.norm("l2")
            resid0 = res

            if(MPI.rank(comm) == 0):
                logger.debug("Iteration: %d, Residual: %.3e, Relative residual: %.3e" % (
                    it, res, rel_res))

            dww = w.copy(deepcopy=True)
            dww.vector()[:] = 0.0

            while (rel_res > rel_tol and res > abs_tol) and it < maxiter:
                it += 1

                A, b = assemble_system(Jac, -Ftotal, bcs)
                solve(A, dww.vector(), b)
                w.vector().axpy(1.0, dww.vector())

                B = assemble(Ftotal)
                for bc in bcs:
                    bc.apply(B)

                rel_res = B.norm("l2")/resid0
                res = B.norm("l2")

                if(MPI.rank(comm) == 0):
                    logger.debug("Iteration: %d, Residual: %.3e, Relative residual: %.3e" % (
                        it, res, rel_res))
                if not(isinstance(res,(float,int))) or not(isinstance(rel_res,(float,int))):
                    logger.critical("Iteration: %d, Residual: %.3e, Relative residual: %.3e" % (
                        it, res, rel_res))
                    raise RuntimeError("Failed Convergence")

            if(rel_res > rel_tol and res > abs_tol):
                raise RuntimeError("Failed Convergence")
