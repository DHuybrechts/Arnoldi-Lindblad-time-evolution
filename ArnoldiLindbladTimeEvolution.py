def arnoldi_time_evolution(H, c_ops, vec_init,n , T, numsteps=100, condition='steady_state', tau=10**(-3), min_check=100, how_often=20):
    """    Computes the orthonormal basis of the (n + 1)-Krylov subspace of 
    mat spanned by {vec, mat*vec, ..., mat^n*vec}.

    Arguments
      H: Hamiltonian of the system
      c_ops: dissipation operators of the systems
      vec_init:  initial density matrix 
      n: maximal dimension of the Krylov subspace, must be an integer
      T: evolution time per iteration
      numsteps: number of steps in the time evolution up to time T
      condition: condition for the wanted eigenvalues; function or keyword
      tau: precision
      min_check: iteration number where checking for convergency starts
      how_often: interval of number of iterations where it checks for convergence
    
    Returns
      Q: m x (n + 1) array of Qobj, the columns are an orthonormal basis of the
        Krylov subspace.
      h: (n + 1) x n array, E on basis S. It is upper Hessenberg and its 
          eigenvalues and eigenvectors converge towards those of E. 
      vals_good: converged eigenvalues of the evolution operator.
      vecs_eff_good: converged right eigenvectors of the evolution operator/
          Liouvillian.
      vals_liouv_good: converged eigenvalues of the Liouvillian.
      
    """
    
    HH=H.full()
    cc_ops=[cc.full() for cc in c_ops]
    def do_liouvillian(rho):
        Lrho= -1.j*(np.dot(HH,rho) - np.dot(rho,HH))
        for jump in cc_ops:
            Lrho=Lrho + np.dot(jump,np.dot(rho,np.conj(np.transpose(jump))))- 0.5*(np.dot(rho,np.dot(np.conj(np.transpose(jump)),jump))+ np.dot(np.dot(np.conj(np.transpose(jump)),jump),rho))
        return Lrho
    
    
    if condition=='steady_state':
        def condition(vv, vals):
            return np.abs(vv)>0.95
    elif isinstance(condition, int):
        NN=condition
        def condition(vv, vals):
            index=np.argsort(abs(vals))
            return vv in vals[index[-NN:]]
        

    q = operator_to_vector(vec_init).full()[:,0]
    q/=np.linalg.norm(q)
    m = q.shape[0]
    sqrtm = int(np.sqrt(m))
    h=np.zeros((n + 1, n), dtype=np.complex128)
    Q=np.zeros((m, n+1), dtype=np.complex128)
    Q[:, 0] = q  # Use it as the first Krylov vector

    tlist=np.linspace(0, T, numsteps)
    
    for k in range(n):
                     
        q=q.reshape((sqrtm,sqrtm))
        q=Qobj(q, dims=vec_init.dims, shape=vec_init.shape)
        sol=mesolve(H, q, tlist, c_ops) #When not using QUTIP, use your own time evolution function here.
                
        v = sol.states[-1]  # Generate a new candidate vector
        v=(v.full()).reshape((m,))
        for j in range(k + 1):  # Subtract the projections on previous vectors
            h[j, k] = np.dot(np.conj(Q[:, j]),v)
            v = v - h[j, k] * Q[:, j]

        h[k + 1, k] = np.linalg.norm(v)
        """
        Orthonormalization is sub-optimal but more numerically stable
        """
        eps = 1e-12  # If v is shorter than this threshold it is the zero vector
        if h[k + 1, k] > eps:  # Add the produced vector to the list, unless
            q = v / h[k + 1, k]  # the zero vector is produced.
            Q[:, k + 1] = q

         
        
        convergence=0
        
        """
        The following part is extremely sub-obptimal. How can we improve it?
        """
        
        if k%how_often==0 and k>min_check:

            vals_eff, vecs_eff=np.linalg.eig(h[0:k, 0:k]) #calculate the effective eigenvalues and eigenvectors

            vals_good=[]
            vals_liouv_good=[]
            vecs_eff_good=[]
            convergence=1
            for jj in range(len(vals_eff)): #loop to see which eigenvalues/ eigenvectors have converged
                if condition(vals_eff[jj], vals_eff)==True:
                    
                    vec=vecs_eff[:,jj]
                    vec=np.dot(Q[:,0:k],vec) #transforms the vector in the effective basis to one expressed in the original basis
                    vec=vec/np.linalg.norm(vec)
                    vec_r=vec.reshape((sqrtm,sqrtm))
                    
                    qq=Qobj(vec_r, dims=vec_init.dims, shape=vec_init.shape)

                    sol=mesolve(H, qq, tlist, c_ops) #When not using QUTIP, use your own time evolution function here.
                    
                    
                    diff=vals_eff[jj]*qq-sol.states[-1]
                    if diff.norm()<tau: #convergence measure
                        convergence*=1
                        vals_good.append(vals_eff[jj]) #appends converged eigenvalue of the evolution operator
                        qq=sol.states[-1].full() 
                        qq/=np.linalg.norm(qq) 
                        Lv=do_liouvillian(qq) 
                        vals_liouv_good.append(np.dot(np.conj(qq).reshape((m,)), Lv.reshape((m,))))   #calculates and appends the eigenvalue of the Liouvillian
                        vecs_eff_good.append(qq) #appends the converged eigenvector
                    else:
                        convergence*=0
                    
            
        if convergence: break
            
        if h[k + 1, k] < eps:  # If that happens, stop iterating.
            break
    print('Total number of iterations = '+str(k))
    return Q, h, vals_good, vecs_eff_good, vals_liouv_good