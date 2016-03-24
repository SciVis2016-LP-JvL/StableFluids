package StableFluids;

import jv.object.PsDebug;

/**
 * Project visualizing smoke like [Jos Stam 1999].
 * 
 * @author		Lukas Polthier, Johannes von Lindheim
 * @version		04.03.2016
 * 				update of FluidSolver, now everything works properly, Gauﬂ-Seidel as linear solver
 */

public class FluidSolver implements Cloneable
{
	boolean colorON = false;
    int n; //width of grid point array
    int m; //height of grid point array
    int size; //number of all grid points
    float dt; //time step

    //constants that determine the nature of the fluid
    //values in parenthesis [0.0f] give the default values of each constant, or the domain
    float visc = 0.0f; //[] constant that determines the viscosity of the fluid
    float diff = 0.0f; //[] constant that determines the effect of diffusion, i.e. how much the density diffuses
    float vorticity = 0.35f; //[0.0f to 1.0f] factor that determines the effect of Vorticity Confinement, i.e. how "curly" the fluid should look like
    						 //the bigger the vorticity, the more curly, but the more crispy
    float weightAir = 5.0f; //[0.025f] constant that determines the upward force due to heat of smoke
    float weightSmoke = 0.4f; //[0.000625f] constant that determines the downward force due to weight of smoke
    float buoyancy = 0.5f;

    float[] d, dOld; //density, old density
    float[] u, uOld; //x-component of the velocity vector field
    float[] v, vOld; //y-component of the velocity vector field
    
    float[] d2, d3, d2Old, d3Old;
    
    float[] tmp;
    float[] curl;
    
    public FluidSolver clone() throws CloneNotSupportedException
    {
    	FluidSolver copy = null;
    	try
    	{
    		copy = (FluidSolver) super.clone();
    		copy.d = d.clone();
    		copy.v = v.clone();
    		copy.u = u.clone();
    		copy.dOld = dOld.clone();
    		copy.vOld = vOld.clone();
    		copy.uOld = uOld.clone();
    		if(colorON)
    		{
	    		copy.d2 = d2.clone();
	    		copy.d2Old = d2Old.clone();
	    		copy.d3 = d3.clone();
	    		copy.d3Old = d3Old.clone();
    		}
    	} catch(CloneNotSupportedException e) {}

    	return copy;
    }
    
    public void setup(int n, int m, float dt, boolean colorON)
    {
    	this.colorON = colorON;
        this.n = n;
        this.m = m;
        this.dt = dt;
        size = (n + 2) * (m + 2); //number of grid points, boundary gets additional row and column
        reset();
    }
    

    public void reset()
    {
        d    = new float[size];
        dOld = new float[size];
        u    = new float[size];
        uOld = new float[size];
        v    = new float[size];
        vOld = new float[size];
        
        if(colorON) {
        	d2    = new float[size];
            d2Old = new float[size];
            d3    = new float[size];
            d3Old = new float[size];
        }
        
        curl = new float[size];
        tmp = new float[size];
                        
    	//initialize velocity and density to 0
        for (int i = 0; i < size; i++)
        {
            u[i] = uOld[i] = v[i] = vOld[i] = 0.0f;
            d[i] = dOld[i] = 0.0f;
        }
        if(colorON) {
        	for (int i = 0; i < size; i++)
            {
        		d2[i] = d2Old[i] = 0.0f;
            	d3[i] = d3Old[i] = 0.0f; 
            }
        }
    }
    
    public void setVisc(float viscosity)
    {
    	this.visc = viscosity;
    }
    
    public float getVisc()
    {
    	return visc;
    }
    
    public void setDiff(float diffusion)
    {
    	this.diff = diffusion;
    }
    
    public float getDiff()
    {
    	return diff;
    }
    
    public void setVorticity(float rotation)
    {
    	this.vorticity = rotation;
    }
    
    public float getVorticity()
    {
    	return vorticity;
    }
    
    public void setWeightAir(float buoyancy)
    {
    	this.weightAir = buoyancy;
    }
    
    public float getWeightAir()
    {
    	return weightAir;
    }
    
    public void setWeightSmoke(float weight)
    {
    	this.weightAir = weight;
    }
    
    public float getWeightSmoke()
    {
    	return weightAir;
    }
    
    public void setBuoyancy(float buoyancy)
    {
    	this.buoyancy = buoyancy;
    }
    
    
    public void resizeArray(int width, int height)
    {
    	int nOld, mOld;
    	size = (n+2)*(m+2);
    	float[] uT = new float[size];
    	float[] vT = new float[size];
    	float[] dT = new float[size];
    	float[] uOldT = new float[size];
    	float[] vOldT = new float[size];
    	float[] dOldT = new float[size];
    	float[] d2OldT = new float[size];
   		float[] d3OldT = new float[size];
   		float[] d2T = new float[size];
   		float[] d3T = new float[size];
    	
    	//store the current value of u, v, d, uOld, vOld, dOld
    	for (int i = 1; i <= Math.min(n, width); i++)
        {
            for (int j = 1; j <= Math.min(m, height); j++)
            {
            	uT[I(i,j)] = u[I(i,j)];
            	vT[I(i,j)] = v[I(i,j)];
            	dT[I(i,j)] = d[I(i,j)];
            	uOldT[I(i,j)] = uOld[I(i,j)];
            	vOldT[I(i,j)] = vOld[I(i,j)];
            	dOldT[I(i,j)] = dOld[I(i,j)];
            }
        }
    	if(colorON) {
    		for (int i = 1; i <= Math.min(n, width); i++)
            {
                for (int j = 1; j <= Math.min(m, height); j++)
                {
                	d2OldT[I(i,j)] = d2Old[I(i,j)];
            		d3OldT[I(i,j)] = d3Old[I(i,j)];
            		d2T[I(i,j)] = d2[I(i,j)];
            		d3T[I(i,j)] = d3[I(i,j)];
                }
            }
    	}
    	//resize the array u,v,d
    	size = (width +2) * (height +2);
    	nOld = n;
    	mOld = m;
    	n = width;
    	m = height;
    	u = new float[size];
    	v = new float[size];
    	d = new float[size];
    	if(colorON) {
    		d2 = new float[size];
    		d3 = new float[size];
    	}
    	uOld = new float[size];
    	vOld = new float[size];
    	dOld = new float[size];
    	if(colorON) {
    		d2Old = new float[size];
    		d3Old = new float[size];
    	}
    	//copy the values back
    	for (int i = 1; i <= Math.min(nOld, n); i++)
        {
            for (int j = 1; j <= Math.min(mOld, m); j++)
            {
            	u[I(i,j)] = uT[i + (nOld + 2) * j];
            	v[I(i,j)] = vT[i + (nOld + 2) * j];
            	d[I(i,j)] = dT[i + (nOld + 2) * j];
            	uOld[I(i,j)] = uOldT[i + (nOld + 2) * j];
            	vOld[I(i,j)] = vOldT[i + (nOld + 2) * j];
            	dOld[I(i,j)] = dOldT[i + (nOld + 2) * j];
            }
        }
    	if(colorON) {
    		for (int i = 1; i <= Math.min(nOld, n); i++)
            {
                for (int j = 1; j <= Math.min(mOld, m); j++)
                {
                	d2[I(i,j)] = d2T[i + (nOld + 2) * j];
            		d3[I(i,j)] = d3T[i + (nOld + 2) * j];
            		d2Old[I(i,j)] = d2OldT[i + (nOld + 2) * j];
            		d3Old[I(i,j)] = d3OldT[i + (nOld + 2) * j];
                }
            }
    	}
    }
    
    public void clearArray()
    {
    	for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= m; j++)
            {
            	u[I(i,j)] = 0;
            	v[I(i,j)] = 0;
            	d[I(i,j)] = 0;
            	uOld[I(i,j)] = 0;
            	vOld[I(i,j)] = 0;
            	dOld[I(i,j)] = 0;
            }
        }
    	if(colorON) {
    		for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= m; j++)
                {
                	d2[I(i,j)] = 0;
            		d3[I(i,j)] = 0;
            		d2Old[I(i,j)] = 0;
            		d3Old[I(i,j)] = 0;
                }
            }
    	}
    }
    
    public void setDensity(int[] red, int[] green, int[] blue)
    {
    	for (int i = 0; i <= n+1; i++)
        {
            for (int j = 0; j <= m+1; j++)
            {
            	d[I(i,j)] = red[I(i,j)];
            	dOld[I(i,j)] = 0;
            }
        }
    	if(colorON) {
    		for (int i = 0; i <= n+1; i++)
            {
                for (int j = 0; j <= m+1; j++)
                {
                	d2[I(i,j)] = green[I(i,j)];
            		d3[I(i,j)] = blue[I(i,j)];
            		d2Old[I(i,j)] = 0;
            		d3Old[I(i,j)] = 0;
                }
            }
    	}
    }
       

    /**
     * This is the main routine. The fluid velocity vector field acts according to the Navier-Stokes equations.
     * The precise numerical approximation is base on a paper by Stam 1999.
     * The 2d vector field is split in 2 arrays u, v, namely the x- and y-component.
     **/

    public void velocitySolver()
    {
        //add velocity that was input by mouse
        inputData(u, uOld);
        inputData(v, vOld);
        //project();

        //add vorticity confinement force
        vorticityConfinementNEW();
              
        //project();

        //add buoyancy force (Auftriebskraft)
        buoyancy();
        
        //buoyancy(vOld);
        //addSource(v, vOld);
        //project();

        
    	diffuse();
    	project();
    	
        advect();
        
        project();
        
        /**
        ///**
        // * This block has been used as a test environment.
        // * 
        //reduce density by a fixed factor, i.e. let the density dissipate
        for (int i = 1; i<=n; i++)
    	{
    		for (int j = 1; j<=m; j++)
    		{
    			d[I(i,j)] = d[I(i,j)] * 0.99f;
    		}
    	}
        
        //add density at specific grid points
        int number = 20;
        for (int i = 1; i<=n/number; i++)
    	{
    		for (int j = 1; j<=m/number; j++)
    		{
    			dOld[I(number*i,number*j)] = dOld[I(number*i,number*j)] +2;
    		}
    	}
    	//*
    	//* This block has been used as a test environment.
    	//*/
        
        // clear all input velocities for next frame
        for (int i = 0; i < size; i++){ uOld[i] = 0; vOld[i] = 0; }
    }
    
    /**
     * Hot smoke goes up, heavy smoke goes down.
     * This method calculates in a simplified version the average temperature and the adds some upward/downward force to the fluid
     */
    
    public void buoyancy()
    {
    	float factor = buoyancy;
        float Tamb = 0;
        //determine average temperature
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= m; j++)
            {
            	if(colorON) {
            		if(d2[I(i, j)] <= 255)
	        		{
	            		Tamb += d2[I(i, j)];
	        		} else {
	        			Tamb = Tamb + 255;
	        		}
            		if(d3[I(i, j)] <= 255)
	        		{
	            		Tamb += d3[I(i, j)];
	        		} else {
	        			Tamb = Tamb + 255;
	        		}
            	}
            	if(d[I(i, j)] <= 255)
        		{
            		Tamb += d[I(i, j)];
        		} else {
        			Tamb = Tamb + 255;
        		}
            }
        }
        Tamb = Tamb / (n * m);

        //hot smoke goes up, heavy smoke goes down
        if(colorON) {
        	for (int i = 1; i <= n; i++)
	        {
	            for (int j = 1; j <= m; j++)
	            {
	                //v[I(i, j)] = v[I(i, j)] + dt * (heat * (d[I(i, j)]) - weight * (d[I(i, j)] - Tamb));
	            	if(d[I(i, j)] + d2[I(i, j)] + d3[I(i, j)] <= 255 + 255 + 255)
	        		{
	            		// v = v + dt * Auftrieb - dt * Gewicht;
	            		v[I(i, j)] = v[I(i, j)] - factor * dt * ( Tamb/255 * weightSmoke + (1-Tamb/255) * weightAir) + factor * dt * ((d[I(i, j)] + d2[I(i, j)] + d3[I(i, j)])/255 * weightSmoke + (1-(d[I(i, j)] + d2[I(i, j)] + d3[I(i, j)])/255) * weightAir);
	        		} else {
	        			v[I(i, j)] = v[I(i, j)] - factor * dt * ( Tamb/255 * weightSmoke + (1-Tamb/255) * weightAir) + factor * dt * (1 * weightSmoke);
	        		}
	            }
	        }
        } else {
	        for (int i = 1; i <= n; i++)
	        {
	            for (int j = 1; j <= m; j++)
	            {
	                //v[I(i, j)] = v[I(i, j)] + dt * (heat * (d[I(i, j)]) - weight * (d[I(i, j)] - Tamb));
	            	if(d[I(i, j)] <= 255)
	        		{
	            		// v = v + dt * Auftrieb - dt * Gewicht;
	            		v[I(i, j)] = v[I(i, j)] - factor * dt * ( Tamb/255 * weightSmoke + (1-Tamb/255) * weightAir) + factor * dt * (d[I(i,j)]/255 * weightSmoke + (1-d[I(i,j)]/255) * weightAir);
	        		} else {
	        			v[I(i, j)] = v[I(i, j)] - factor * dt * ( Tamb/255 * weightSmoke + (1-Tamb/255) * weightAir) + factor * dt * (1 * weightSmoke);
	        		}
	            }
	        }
        }
    }
    
    
    /**
     * The vorticity confinement step adds some small-scale details to the flow that have been lost by the coarse grid.
     * First we compute the rotation or curl of the 2d-vector field. The curl, being defined as $\nabla x u$ is a 3d vector field
     * that is perpendicular to the xy-plane. Thus we refer to the curl as a 1d array.
     * Then we calculate the normalized gradient of the absolute curl, i.e. the normal vector field that points towards points with high vorticity.
     * Finally we add the force given by $\epsilon h (N x  curl)$ to our original vector field, where N is the normalized gradient of the absolute curl.
     * I.e. we add some force that is perpendicular to the gradient of the absolute curl and scaled with the curl.
     */
    
    public void vorticityConfinementNEW()
    {
    	
    	float n1,n2,factor;
    	float[] rot;
    	rot = new float[size];
    	//first compute the absolute curl at each point of the grid
    	for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= m; j++)
            {
                rot[I(i, j)] = Math.abs( (u[I(i, j + 1)] - u[I(i, j - 1)]) - (v[I(i + 1, j)] - v[I(i - 1, j)]) ) /2;
            }
        }
    	
    	//then compute the resulting force
    	for (int i = 2; i < n; i++)
        {
            for (int j = 2; j < m; j++)
            {
            	//compute the gradient of the absolute curl
            	//(note that we drop the factor n as we normalize anyway)
            	n1 = (rot[I(i+1,j)] - rot[I(i-1,j)] ) / 2 ;
            	n2 = (rot[I(i,j+1)] - rot[I(i,j-1)] ) / 2 ;
            	
            	//normalize the gradient
            	factor = (float) Math.sqrt(n1 * n1 + n2 * n2) + 0.000000001f; //length
            	n1 = n1 / factor;
            	n2 = n2 / factor;
            	
            	//add force to the vector field
            	u[I(i,j)] = u[I(i,j)] - n2 * vorticity  * ((u[I(i, j + 1)] - u[I(i, j - 1)]) - (v[I(i + 1, j)] - v[I(i - 1, j)]))/2;
            	v[I(i,j)] = v[I(i,j)] + n1 * vorticity  * ((u[I(i, j + 1)] - u[I(i, j - 1)]) - (v[I(i + 1, j)] - v[I(i - 1, j)]))/2;
            }
        }
    }
    
    
    /**
     * Do the advection step, i.e. let the ether "float" a little, by using the particle tracer P : \Omega \times \mathbb{R} \rightarrow \Omega.
     * Starting at $x\in \Omega$, $P(x,-\delta t)$ gives the position at previous time $-\delta t$ by going backwards through the velocity field.
     * 
     * In practice, for all grid cells, we trace the grid cells center backwards through the velocity field.
     * Then we interpolate from the grid of the previous timestep the corresponding vector at this position.
     * We assign this value to the current grid cell.
     * 
     * For simplicity we handle the x-and y-component of the velocity field separately.
	 *
     **/
    
    private void advect()
    {
		float x, y; //factors in the convex combination
		int xI, yI;
		for (int i = 1; i<=n; i++)
    	{
    		for (int j = 1; j<=m; j++)
    		{
    			uOld[I(i,j)] = u[I(i,j)];
    			vOld[I(i,j)] = v[I(i,j)];
    		}
    	}
		
    	for (int i = 1; i<=n; i++)
    	{
    		for (int j = 1; j<=m; j++)
    		{
    			//go backwards through the velocity field
    			x = i - uOld[I(i,j)] * dt * m;
    			y = j - vOld[I(i,j)] * dt * m;
    			
    			//interpolate the corresponding velocity vector at the position P(x,-\delta t)
    			//catch boundary overshooting
    			if (x > n + 0.5) x = n + 0.5f;
                if (x < 0.5)     x = 0.5f;
                if (y > m + 0.5) y = m + 0.5f;
                if (y < 0.5)     y = 0.5f;
                //compute convex-combination to get the weights                
    			xI = (int) x;
    			yI = (int) y;
    			
    			//assign the resulting vector to the output array
    			u[I(i,j)] = (1-x+xI) * ( (1-y+yI) * uOld[I(xI,yI)] + (y-yI) * uOld[I(xI,yI+1)] )
    					+ (x-xI) * ( (1-y+yI) * uOld[I(xI+1,yI)] + (y-yI) * uOld[I(xI+1,yI+1)] );
    			v[I(i,j)] = (1-x+xI) * ( (1-y+yI) * vOld[I(xI,yI)] + (y-yI) * vOld[I(xI,yI+1)] )
    					+ (x-xI) * ( (1-y+yI) * vOld[I(xI+1,yI)] + (y-yI) * vOld[I(xI+1,yI+1)] );
    		}
    	}
    	
    	boundaryCondition(1, u);
    	boundaryCondition(2, v);
    	
    	
    }
    
    
    /**
     * Do the diffusion step, i.e. we consider the friction of the ether with itself.
     * Thus we need to solve the implicit system $(Id - dt * \nu \delta)w_3 = w_2$.
     * This leads to the problem of solving a sparse linear system.
     * Method to solve the sparse linear system of choice is the Gauss-Seidel method.
     */
    
    private void diffuse()
    {
		float factor = visc * dt * m * m;
		// PsDebug.message("factor:"+ Float.toString(factor));
		float[] neu;
		float[] neu2;
		neu = new float[size];
		neu2 = new float[size];
		neu = uOld;
		neu2 = vOld;
		
    	for (int k = 1;k< 20; k++)
    	{
    		for (int i = 1;i<=n; i++)
    		{
    			for (int j = 1;j<=m; j++)
    			{
    				neu[I(i,j)] = ( u[I(i,j)] +factor*neu[I(i+1,j)] +factor*neu[I(i,j+1)] +factor*neu[I(i-1,j)] +factor*neu[I(i,j-1)] ) / ( 1 + 4 * factor);
    				neu2[I(i,j)] = ( v[I(i,j)] +factor*neu2[I(i+1,j)] +factor*neu2[I(i,j+1)] +factor*neu2[I(i-1,j)] +factor*neu2[I(i,j-1)] ) / ( 1 + 4 * factor);
    			}
    		}    		
    	}    
    	for (int i = 1;i<=n; i++)
		{
			for (int j = 1;j<=m; j++)
			{
				u[I(i,j)] = neu[I(i,j)];
				v[I(i,j)] = neu2[I(i,j)];
			}
		}
    }
    
    
    /**
     * Do the projection step. This ensures that the resulting velocity vector field is divergence-free, i.e. mass preserving.
     * As in the diffusion step we need to solve a sparse linear system:
     * By Hodge, $w_3 = u + \nabla q$ where $div u = 0$ and $q$ is a scalar field.
     * The decomposition is implicitly given by the solution of the Poisson equation $\delta q = div w_3$.
     * First we solve the Poisson equation with the difference method and Gauﬂ-Seidel iteration,
     * then we compute the gradient field $\nabla q$.
     * Finally the new vector field $w_4 = w_3 - \nabla q$. 
     */
    
    void project()
    {
		float[] neu2;
		float[] div;
		neu2 = new float[size];
		div = new float[size];
		//compute the divergence of the vector field
    	for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= m; j++)
            {
            	div[I(i,j)] = ( u[I(i+1,j)] - u[I(i-1,j)] ) / (2 * n) + ( v[I(i,j+1)] - v[I(i,j-1)] ) / (2 * n);
            	neu2[I(i,j)] = 0;
            }
        }
    	boundaryCondition(0, div);
    	boundaryCondition(0, neu2);
    	
    	
    	//solve the resulting Poisson equation
    	for (int k = 1; k < 20; k++)
    	{
	    	for (int i = 1; i <= n; i++)
	        {
	            for (int j = 1; j <= m; j++)
	            {
	            	neu2[I(i,j)] = ( div[I(i,j)] - neu2[I(i-1,j)] - neu2[I(i,j-1)] - neu2[I(i+1,j)] - neu2[I(i,j+1)]) / (-4); 
	            }
	        }
    	}
    	boundaryCondition(0, neu2);

    	//finally we obtain the mass preserving field by subtracting the gradient of neu from the original vector field
    	for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= m; j++)
            {
            	u[I(i,j)] = u[I(i,j)] - (neu2[I(i+1,j)] - neu2[I(i-1,j)]) /2 * n;
            	v[I(i,j)] = v[I(i,j)] - (neu2[I(i,j+1)] - neu2[I(i,j-1)]) /2 * n;
            }
        }
    	boundaryCondition(1, u);
        boundaryCondition(2, v);
    }

    
    public void densitySolver()
    {
    	float[] neu;
		neu = new float[size];
        // add density inputed by mouse
        inputData(d, dOld);
        if(colorON) {
        	inputData(d2, d2Old);
        	inputData(d3, d3Old);
        }
        
        //let the density diffuse
        float factor = diff * dt * m * m ;
        neu = dOld;
        for (int k = 1;k<= 20; k++)
    	{
    		for (int i = 1;i<=n; i++)
    		{
    			for (int j = 1;j<=m; j++)
    			{
    				neu[I(i,j)] = ( d[I(i,j)] +factor*neu[I(i+1,j)] +factor*neu[I(i,j+1)] +factor*neu[I(i-1,j)] +factor*neu[I(i,j-1)] ) / ( 1 + 4 * factor);
    			}
    		}    		
    	} 
        dOld = neu;
        
        if(colorON) {
    		neu = new float[size];
            
            //let the density diffuse
            neu = d2Old;
            for (int k = 1;k<= 20; k++)
        	{
        		for (int i = 1;i<=n; i++)
        		{
        			for (int j = 1;j<=m; j++)
        			{
        				neu[I(i,j)] = ( d2[I(i,j)] +factor*neu[I(i+1,j)] +factor*neu[I(i,j+1)] +factor*neu[I(i-1,j)] +factor*neu[I(i,j-1)] ) / ( 1 + 4 * factor);
        			}
        		}    		
        	} 
            d2Old = neu;
            
    		neu = new float[size];
            
            //let the density diffuse
            neu = d3Old;
            for (int k = 1;k<= 20; k++)
        	{
        		for (int i = 1;i<=n; i++)
        		{
        			for (int j = 1;j<=m; j++)
        			{
        				neu[I(i,j)] = ( d3[I(i,j)] +factor*neu[I(i+1,j)] +factor*neu[I(i,j+1)] +factor*neu[I(i-1,j)] +factor*neu[I(i,j-1)] ) / ( 1 + 4 * factor);
        			}
        		}    		
        	} 
            d3Old = neu;
        }
        
        //let the density advect
        float x, y, f1, f2, f3, f4 = 0;
		int xI, yI = 0;
    	for (int i = 1; i<=n; i++)
    	{
    		for (int j = 1; j<=m; j++)
    		{
    			//go backwards through the velocity field
    			x = i - u[I(i,j)] * dt * m;
    			y = j - v[I(i,j)] * dt * m;
    			
    			//interpolate the corresponding velocity vector at the position P(x,-\delta t)
    			//catch boundary overshooting
    			if (x > n + 0.5) x = n + 0.5f;
                if (x < 0.5)     x = 0.5f;
                if (y > m + 0.5) y = m + 0.5f;
                if (y < 0.5)     y = 0.5f;
                //compute convex-combination to get the weights                
    			xI = (int) x;
    			yI = (int) y;
    			
    			f1 = x - xI;
    			f2 = 1 - f1;
    			f3 = y - yI;
    			f4 = 1 - f3;
    			
    			//assign the resulting vector to the output array
    			d[I(i,j)] = f2 * ( f4 * dOld[I(xI,yI)] + f3 * dOld[I(xI,yI+1)])
    				  + f1 * ( f4 * dOld[I(xI+1,yI)] + f3 * dOld[I(xI+1,yI+1)]);
    		}
    		boundaryCondition(0, d);
    	}      
    	
    	if(colorON) {
    		for (int i = 1; i<=n; i++)
        	{
        		for (int j = 1; j<=m; j++)
        		{
        			//go backwards through the velocity field
        			x = i - u[I(i,j)] * dt * m;
        			y = j - v[I(i,j)] * dt * m;
        			
        			//interpolate the corresponding velocity vector at the position P(x,-\delta t)
        			//catch boundary overshooting
        			if (x > n + 0.5) x = n + 0.5f;
                    if (x < 0.5)     x = 0.5f;
                    if (y > m + 0.5) y = m + 0.5f;
                    if (y < 0.5)     y = 0.5f;
                    //compute convex-combination to get the weights                
        			xI = (int) x;
        			yI = (int) y;
        			
        			f1 = x - xI;
        			f2 = 1 - f1;
        			f3 = y - yI;
        			f4 = 1 - f3;
        			
        			//assign the resulting vector to the output array
        			d2[I(i,j)] = f2 * ( f4 * d2Old[I(xI,yI)] + f3 * d2Old[I(xI,yI+1)])
        				  + f1 * ( f4 * d2Old[I(xI+1,yI)] + f3 * d2Old[I(xI+1,yI+1)]);
        		}
        		boundaryCondition(0, d2);
        	}
    		for (int i = 1; i<=n; i++)
        	{
        		for (int j = 1; j<=m; j++)
        		{
        			//go backwards through the velocity field
        			x = i - u[I(i,j)] * dt * m;
        			y = j - v[I(i,j)] * dt * m;
        			
        			//interpolate the corresponding velocity vector at the position P(x,-\delta t)
        			//catch boundary overshooting
        			if (x > n + 0.5) x = n + 0.5f;
                    if (x < 0.5)     x = 0.5f;
                    if (y > m + 0.5) y = m + 0.5f;
                    if (y < 0.5)     y = 0.5f;
                    //compute convex-combination to get the weights                
        			xI = (int) x;
        			yI = (int) y;
        			
        			f1 = x - xI;
        			f2 = 1 - f1;
        			f3 = y - yI;
        			f4 = 1 - f3;
        			
        			//assign the resulting vector to the output array
        			d3[I(i,j)] = f2 * ( f4 * d3Old[I(xI,yI)] + f3 * d3Old[I(xI,yI+1)])
        				  + f1 * ( f4 * d3Old[I(xI+1,yI)] + f3 * d3Old[I(xI+1,yI+1)]);
        		}
        		boundaryCondition(0, d3);
        	}
    	}

        // clear input density array for next frame
        for (int i = 0; i < size; i++) dOld[i] = 0;
        if(colorON) {
        	for (int i = 0; i < size; i++) d2Old[i] = 0;
        	for (int i = 0; i < size; i++) d3Old[i] = 0;
        }
    }

    
    /**
     * This method simply adds the second array to the first array.
     * Used to add 
     * @param x the original array, to which we add the input data
     * @param x0 the additional force / input data
     */
    
    private void inputData(float[] x, float[] x0)
    {
        for (int i = 0; i < size; i++)
        {
            x[i] += dt * x0[i];
        }
    }

    
    /**
     * This function handles the boundary conditions.
     * @param b flag that gives the behavior at the boundary
     * 		b = 0 -> no reflection at any boundary
     * 		b = 1 -> reflection on the left and right boundary
     * 		b = 2 -> reflection on the upper and lower boundary
     * @param x array that should interact with the boundary
     */
    
    private void boundaryCondition(int b, float[] x)
    {
    	//go along the boundary points
        for (int i = 1; i <= n; i++)
        {          
            //upper and lower boundary
            x[I(  i, 0  )] = b == 2 ? -x[I(i, 1)] : x[I(i, 1)];
            x[I(  i, m+1)] = b == 2 ? -x[I(i, m)] : x[I(i, m)];
        }
        
        for (int i = 1; i <= m; i++)
        {
    		//boundary on the left and on the right
            x[I(  0, i  )] = b == 1 ? -x[I(1, i)] : x[I(1, i)];
            x[I(n+1, i  )] = b == 1 ? -x[I(n, i)] : x[I(n, i)];
        }

    	//the value at corners is the average of the effect of both boundaries
        try {
        x[I(  0,   0)] = 0.5f * (x[I(1, 0  )] + x[I(  0, 1)]);
        }
        catch (NullPointerException e)
        {
        	// System.out.println("nullpointerexception! x.length == " + String.valueOf(x.length));
        }
        x[I(  0, m+1)] = 0.5f * (x[I(1, m+1)] + x[I(  0, m)]);
        x[I(n+1,   0)] = 0.5f * (x[I(n, 0  )] + x[I(n+1, 1)]);
        x[I(n+1, m+1)] = 0.5f * (x[I(n, m+1)] + x[I(n+1, m)]);

    }

    // method for indexing 1d arrays
    private int I(int i, int j){ return i + (n + 2) * j; }
}