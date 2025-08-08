from dolfin import *
import logging
logger = logging.getLogger(__name__)
import math

class NSolver(object):

    def __init__(self, params):

        self.parameters = params
        self.isfirstiteration = 0
        self.default_parameters={"rel_tol": 1e-6,
                                "abs_tol": 1e-6,
                                "max_iter": 5}
        for key in self.default_parameters:
            if key in self.parameters:
                self.default_parameters[key]=self.parameters[key]

    def solvenonlinear(self):
        if "lr" in self.parameters:
            lr = self.parameters["lr"]
        else:
            lr = 1.0
        logger.debug("Learning rate for solver set at "+repr(lr))
        abs_tol = self.default_parameters["abs_tol"]
        rel_tol = self.default_parameters["rel_tol"]
        maxiter = self.default_parameters["max_iter"]
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
            
            for FF in ["1","2","3_Atr","3_Ven","4","5","6"]:
                if 'F'+FF in self.parameters:
                    if self.parameters['F'+FF] is not None:
                        Ftotaltemp = assemble(self.parameters['F'+FF])
                        if isinstance(Ftotaltemp,float):
                            logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp))
                        else:
                            logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp.norm("l2")))
            
            prev_res=float('inf')
            saved_w_vector = w.copy(deepcopy=True)
            if(self.isfirstiteration == 0):
                B = assemble(Ftotal)
                for bc in bcs:
                    bc.apply(B)
                rel_res = 1.0
                res = B.norm("l2")
                dww = w.copy(deepcopy=True)
                dww.vector()[:] = 0.0
                if(MPI.rank(comm) == 0):
                    logger.debug("Start. Residual: %.3e, Relative residual: %.3e" % (
                        res, rel_res))
                A, b = assemble_system(Jac, -Ftotal, bcs)
                resid0 = b.norm("l2")
                rel_res = b.norm("l2")/resid0
                res = resid0
                prev_res=res
                if(MPI.rank(comm) == 0):
                    logger.debug("isfirstiteration")
                    logger.debug("Iteration: %d, Residual: %.3e, Relative residual: %.3e" % (
                        it, res, rel_res))
                solve(A, dww.vector(), b)
                w.vector().axpy(lr, dww.vector())
                for debug_call in self.parameters["debug_call"]:
                    debug_call[0](*debug_call[1:])
                it += 1

            self.isfirstiteration = 1

            B = assemble(Ftotal)
            for bc in bcs:
                bc.apply(B)

            rel_res = 1.0
            res = B.norm("l2")
            resid0 = res
            while res>(2.*prev_res) or math.isnan(res):
                if math.isnan(res):
                    for FF in ["1","2","3_Atr","3_Ven","4","5","6"]:
                        if 'F'+FF in self.parameters:
                            if self.parameters['F'+FF] is not None:
                                Ftotaltemp = assemble(self.parameters['F'+FF])
                                if isinstance(Ftotaltemp,float):
                                    logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp))
                                else:
                                    logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp.norm("l2")))
                w.vector()[:]=saved_w_vector.vector()[:]
                #w.vector().axpy(-lr, dww.vector())
                logger.debug("reverted back")
                B = assemble(Ftotal)
                for bc in bcs:
                    bc.apply(B)

                rel_res = 1.0
                res = B.norm("l2")
                resid0 = res
                for FF in ["1","2","3_Atr","3_Ven","4","5","6"]:
                    if 'F'+FF in self.parameters:
                        if self.parameters['F'+FF] is not None:
                            Ftotaltemp = assemble(self.parameters['F'+FF])
                            if isinstance(Ftotaltemp,float):
                                logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp))
                            else:
                                logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp.norm("l2")))
                
                lr=lr*0.9
                if(MPI.rank(comm) == 0):
                    logger.debug("Reduced learning rate to: %.3e" % (lr))
                if lr<1e-8:
                    raise RuntimeError("Failed Convergence: learning rate too low.")
                w.vector().axpy(lr, dww.vector())
                for debug_call in self.parameters["debug_call"]:
                    debug_call[0](*debug_call[1:])
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
            saved_w_vector = w.copy(deepcopy=True)
            dww.vector()[:] = 0.0
            prev_res=res
            reduce_lr_condition=2.

            while (rel_res > rel_tol and res > abs_tol) and it < maxiter:
                '''
                for FF in ["1","2","3_Atr","3_Ven","4","5","6"]:
                    if 'F'+FF in self.parameters:
                        if self.parameters['F'+FF] is not None:
                            Ftotaltemp = assemble(self.parameters['F'+FF])
                            if isinstance(Ftotaltemp,float):
                                logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp))
                            else:
                                logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp.norm("l2")))
                '''
                it += 1

                A, b = assemble_system(Jac, -Ftotal, bcs)
                solve(A, dww.vector(), b)
                w.vector().axpy(lr, dww.vector())
                for debug_call in self.parameters["debug_call"]:
                    debug_call[0](*debug_call[1:])
                B = assemble(Ftotal)
                for bc in bcs:
                    bc.apply(B)

                rel_res = B.norm("l2")/resid0
                res = B.norm("l2")
                while res>(reduce_lr_condition*prev_res)  or math.isnan(res):
                    if math.isnan(res):
                        for FF in ["1","2","3_Atr","3_Ven","4","5","6"]:
                            if 'F'+FF in self.parameters:
                                if self.parameters['F'+FF] is not None:
                                    Ftotaltemp = assemble(self.parameters['F'+FF])
                                    if isinstance(Ftotaltemp,float):
                                        logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp))
                                    else:
                                        logger.debug(FF+" : Residual: %.3e" % (Ftotaltemp.norm("l2")))
                    w.vector()[:]=saved_w_vector.vector()[:]
                    #w.vector().axpy(-lr, dww.vector())
                    lr=lr*0.9
                    if(MPI.rank(comm) == 0):
                        logger.debug("Reduced learning rate to: %.3e" % (lr))
                    if lr<1e-8:
                        logger.critical("Failed Convergence: learning rate too low.")
                        break
                    w.vector().axpy(lr, dww.vector())
                    for debug_call in self.parameters["debug_call"]:
                        debug_call[0](*debug_call[1:])
                    B = assemble(Ftotal)
                    for bc in bcs:
                        bc.apply(B)

                    rel_res = B.norm("l2")/resid0
                    res = B.norm("l2")
                prev_res=res
                saved_w_vector = w.copy(deepcopy=True)

                if(MPI.rank(comm) == 0):
                    logger.debug("Iteration: %d, Residual: %.3e, Relative residual: %.3e" % (
                        it, res, rel_res))
                
                if not(isinstance(res,(float,int))) or not(isinstance(rel_res,(float,int))) or math.isnan(res) or math.isnan(rel_res):
                    logger.critical("Iteration: %d, Residual: %.3e, Relative residual: %.3e" % (
                        it, res, rel_res))
                    raise RuntimeError("Failed Convergence")
                if it == maxiter and reduce_lr_condition>1.:
                    it=1
                    reduce_lr_condition-=1.
                if lr<1e-6:
                    break
            if(rel_res > rel_tol and res > abs_tol):
                logger.critical("Maximum iteration reached without conergence. Residual: %.3e, Relative residual: %.3e" % (res, rel_res))
