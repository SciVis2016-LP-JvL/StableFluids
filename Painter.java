package StableFluids;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Image;
import java.awt.Point;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.image.MemoryImageSource;

import jv.anim.PsAnimation;
import jv.anim.PsTimeEvent;
import jv.number.PdColor;
import jv.number.PuInteger;
import jv.number.PuDouble;
import jv.object.PsDebug;
import jv.project.PjProject;
import jv.project.PvCameraIf;
import jv.project.PvCameraListenerIf;
import jv.project.PvDisplayIf;
import jv.project.PvPickEvent;
import jv.project.PvCameraEvent;
import jv.project.PvViewerIf;
import jv.vecmath.PdVector;
import jv.vecmath.PiVector;
import jv.objectGui.PsImage;

import jvx.curve.PgBezierCurve;

/**
 * Project visualizing smoke like [Jos Stam 1999].
 * 
 * @author		Lukas Polthier, Johannes von Lindheim, Konrad Polthier (Julia Set project)
 * @version		0.0, 20.12.15 reused JuliaSets project as template for Stable Fluids application
 */
@SuppressWarnings("serial")
public class Painter extends PjProject implements ComponentListener {
	boolean colorON = false;
	boolean colorONset = false;
	boolean colorChange = false;
	boolean isFrozen = false;
	int whichColor = 1;
	// Display of main window
	protected	PvDisplayIf			m_disp;
	// Image to be used as background in display
	protected	Image				m_image;
	// Height of display canvas
	protected	int					m_imageHeight;
	// Save imageHeight
	protected	int					m_oldImageHeight;
	// Width of display canvas
	protected	int					m_imageWidth;
	// Save imageWidth
	protected	int					m_oldImageWidth;
	// BlockWidth, i.e. number of blocks horizontally
	protected	int					m_numBlocksX;
	// BlockHeight, i.e. number of blocks vertically
	protected	int					m_numBlocksY;
	// Producer of image m_image from pixel array
	private		MemoryImageSource	m_mis;
	// Pixel array which stores the image of a textured element. Only used when texture enabled
	private		PiVector			m_pix;
	// Pixel array which stores known used iterations
	private		PdVector			m_density;
	private		PdVector			m_density2;
	private		PdVector			m_density3;

	// Handling mouse input
	private		PiVector			m_forceTraceX;
	private		PiVector			m_forceTraceY;
	protected	PuInteger			m_forceRadius;
	protected	PuDouble			m_forceConst;
	protected	PuDouble			m_buoyancy;
	protected	PuDouble			m_diffusion;
	protected	PuDouble			m_viscosity;
	protected	PuDouble			m_vorticity;
	private		PiVector			m_densityTraceX;
	private		PiVector			m_densityTraceY;
	private		PdVector			m_oldInputDensity;
	protected	PuInteger			m_densityRadius;
	protected	PuDouble			m_densityConst;
	private final int				NUM_PARTS_DENSITY_TRACE_SPLINES = 20;

	// Determines size of uniformly colored pixels blocks in images.
	// Using m_blockSize == 1 leads to highest resolution where each image pixel is really computed.
	protected	PuInteger			m_blockSize;
	protected	PuInteger			m_oldBlockSize;
	
	// Determines, whether computeImage() shall be invoked in setTime
	private		boolean				m_readyForComputeImage;
	private		boolean				m_currentlyResizing;
	
	private		PuDouble			m_time;
	
	// Calculates the simulation data
	private 	FluidSolver			m_fluidSolver;
	private 	FluidSolver			m_oldFluidSolver;

    private 	PuDouble 			m_dt;
    
	public Painter() {
		super("Stable Fluids");
		myReset();
	}
	public void myReset()
	{
		// PsDebug.message("reset");
		m_pix					= new PiVector();
		m_density				= new PdVector();
		m_density2				= new PdVector();
		m_density3				= new PdVector();
		m_forceTraceX			= new PiVector();
		m_forceTraceY			= new PiVector();
		m_forceRadius			= new PuInteger("Force radius");
		m_forceConst			= new PuDouble("Force constant");
		m_buoyancy				= new PuDouble("Buoyancy force", this);
		m_diffusion				= new PuDouble("Diffusion", this);
		m_viscosity				= new PuDouble("Viscosity", this);
		m_vorticity				= new PuDouble("Vorticity", this);
		m_densityTraceX			= new PiVector();
		m_densityTraceY			= new PiVector();
		m_oldInputDensity		= new PdVector();
		m_densityRadius			= new PuInteger("Density radius");
		m_densityConst			= new PuDouble("Density constant");
		m_fluidSolver			= new FluidSolver();
		m_oldFluidSolver		= new FluidSolver();
		m_currentlyResizing		= false;
		m_readyForComputeImage	= true;
		m_dt					= new PuDouble("Timestep");
		m_time					= new PuDouble("Time");
		m_blockSize				= new PuInteger("Block Size", this);
		m_oldBlockSize			= new PuInteger("Old block size", this);
		
		// Initialize animation
		initAnimation();

		if (getClass() == Painter.class) {
			init();
		}
	}

	// Called, when reset-Button is pressed
	public synchronized void buttonReset()
	{
		init();
		initAnimation();
		myStart();
		initFluidSolver();
	}
	
	public void buttonClear()
	{
		m_fluidSolver.clearArray();
	}
	
	public void buttonFlipColor()
	{
		if(colorON) {
			if(whichColor <=5) {
				whichColor = whichColor + 1;
			} else {
				whichColor = 1;
			}
		}
	}
	
	public void buttonImportImage()
	{
		if(colorON) {} else { buttonColorOnOff(); }
		PsImage bild;
		float factor = 0.00392156863f;
		float alpha;
		bild = new PsImage("mysource/StableFluids/test.png");
		// bild = new PsImage("myProjects/StableFluids/test.png");
		//bild.setSize(m_numBlocksX, m_numBlocksY);
		int[] pixelBild;
		Image bild2 = bild.getImage();
		pixelBild = new int[(bild.getHeight()) * (bild.getWidth())];
		pixelBild = PsImage.getPixels( bild2 );
		int[] red, green, blue;
		red = new int[(m_numBlocksX+2) * (m_numBlocksY+2)];
		green = new int[(m_numBlocksX+2) * (m_numBlocksY+2)];
		blue = new int[(m_numBlocksX+2) * (m_numBlocksY+2)];
		
		for(int i=0;i<m_numBlocksX; i++) {
			for(int j=0;j<m_numBlocksY; j++) {
				Color farbe = new Color( pixelBild[ i + bild.getWidth() * (j) ], false );
				//alpha = Math.round( farbe.getAlpha() ) / 255;
				alpha = 1;
				/**
				red[Id(i,j)] = Math.round( (255 - farbe.getRed())*factor / alpha );
				green[Id(i,j)] = Math.round( (255 - farbe.getGreen())*factor / alpha );
				blue[Id(i,j)] = Math.round( (255 - farbe.getBlue())*factor/ alpha );
				*/
				
				
				red[Id(i,j)] = 255 - farbe.getRed();
				green[Id(i,j)] = 255 - farbe.getGreen();
				blue[Id(i,j)] = 255 - farbe.getBlue();
				
				
				/**
				red[Id(i,j)] = 255;
				green[Id(i,j)] = 200;
				blue[Id(i,j)] = 0;
				*/
			}
		}
		//set density in the fluidsolver
		m_fluidSolver.setDensity(red,green,blue);	
	}
	
	public void buttonFreeze()
	{
		if(isFrozen) {
			isFrozen = false;
		} else {
			isFrozen = true;
		}
	}
	
	public void buttonColorOnOff()
	{
		colorChange = true;
		if(colorONset) {
			colorONset = false;
		} else {
			colorONset = true;
		}
		buttonReset();
	}

	public void init()
	{
		super.init();
		colorON = colorONset;

		m_image					= null;
		m_pix.setSize(0);
		m_density.setSize(0);
		m_density2.setSize(0);
		m_density3.setSize(0);
		m_densityRadius.setBounds(1, 40, 1, 5);
		m_densityRadius.setValue(15);
		m_densityConst.setBounds(0.0, 1000.0, 10.0, 100.0);
		m_densityConst.setValue(500.0);
		m_forceTraceX.setSize(0);
		m_forceTraceY.setSize(0);
		m_forceRadius.setBounds(1, 100, 1, 5);
		m_forceRadius.setValue(40);
		m_forceConst.setBounds(0.0, 1.0, 0.01, 0.1);
		m_forceConst.setValue(0.4);
		m_buoyancy.setBounds(-0.5, 0.5, 0.01, 0.1);
		m_buoyancy.setValue(0.1);
		//m_buoyancy.setValue(0.0);
		m_diffusion.setBounds(0.0, 1.0, 0.01, 0.1);
		m_diffusion.setValue(0.0);
		m_viscosity.setBounds(0.0, 1.0, 0.01, 0.1);
		m_viscosity.setValue(0.0);
		m_vorticity.setBounds(0.0, 1.0, 0.01, 0.1);
		m_vorticity.setValue(0.35);
		m_oldInputDensity.setSize(0);
		m_densityTraceX.setSize(0);
		m_densityTraceY.setSize(0);
		resetMouseTraces();
		m_time.setValue(0);
		m_dt.setValue(0.2f);
		
		m_imageHeight			= 0;
		m_imageWidth			= 0;
		m_oldImageHeight		= 0;
		m_oldImageWidth			= 0;
		m_numBlocksX			= 0;
		m_numBlocksY			= 0;
		
		m_blockSize.init();
		m_blockSize.setBounds(1, 10, 1, 2);
		m_blockSize.setValue(1);
		m_oldBlockSize.init();
		m_oldBlockSize.setBounds(1, 10, 1, 2);
		m_oldBlockSize.setValue(1);
	}

	/**
	 * Called when project is launched by viewer on applet start.
	 */
	public void start()
	{
		myStart();
		super.start();
	}
	
	@SuppressWarnings("deprecation")
	private void myStart()
	{
		if (m_disp == null)
			m_disp = getDisp();
		m_disp.showScenegraph(false); // :P
		
		// Adjust sizes of images to dimension of display canvas
		if (resizeImage(m_disp))
		{
			computeImage();
			m_disp.update(null);
		}

		if (m_disp != null)
		{
			m_disp.selectCamera(PvCameraIf.CAMERA_ORTHO_XY);
			m_disp.setBackgroundImageFit(PvDisplayIf.IMAGE_RESIZE);
			m_disp.setMajorMode(PvDisplayIf.MODE_INITIAL_PICK);
			m_disp.showBackgroundImage(true);
		}
		
		// Set time inverval and start the animation
		startAnimation();

		update(this);

	}
	
	/**
	 * Update the class whenever a child has changed.
	 * Method is usually invoked from the children.
	 */
	public boolean update(Object event) 
	{
		// Null
		if (m_disp == null)
			return super.update(event);
		// This
		if (event == this) 
		{
			return super.update(this);
		}
		// BlockSize
		else if (event == m_blockSize) 
		{
			// Update resolution-values and fluidSolver to new blockSize
			changeBlockSize();

			computeImage();
			// m_disp.update(null);
			
			return true;
		}
		// Vorticity confinement
		else if (event == m_vorticity)
		{
			m_fluidSolver.setVorticity((float) m_vorticity.getValue());
			return true;
		}
		// Viscosity
		else if (event == m_viscosity)
		{
			m_fluidSolver.setVisc((float) m_viscosity.getValue() / 8);
			return true;
		}
		// Diffusion
		else if (event == m_diffusion)
		{
			m_fluidSolver.setDiff((float) m_diffusion.getValue() / 100000);
			return true;
		}
		// Buoyancy
		else if (event == m_buoyancy)
		{
			m_fluidSolver.setBuoyancy(10.0f * (float) m_buoyancy.getValue());
			return true;
		}
		return super.update(event);
	}
	
	// Initialize animation: Create time listener, reset time intervals etc.
	private void initAnimation()
	{
		// PsDebug.message("initAnimation");
		if (!super.hasAnimation())
		{
			PsAnimation anim = new PsAnimation();
			anim.setName("Animation");	// Set title of the animation dialog.
			if (!anim.hasTimeListener(this))
			{
				// PsDebug.message("create new time listener");
				anim.addTimeListener(this);
			}
		}

		PsAnimation anim = super.getAnimation();
		m_readyForComputeImage = false;
		anim.setTime(0.0);
		m_readyForComputeImage = true;
		anim.setTimeInterval(0., .2, .1, .1);
	}
	
	// Start animation, if there is one
	private void startAnimation()
	{
		// PsDebug.message("start Animation");
		if (super.hasAnimation())
			super.getAnimation().start();
	}
	
	public boolean setTime(PsTimeEvent timeEvent)
	{
		// PsDebug.message("setTime");
		m_time.setValue(Math.round(timeEvent.getTime()));
		if (m_readyForComputeImage)
			computeImage();
		m_disp.update(null);
		if (super.hasAnimation())
			super.getAnimation().setTimeInterval(0., super.getAnimation().getMaxTime() + 1.);
		return super.update(this);				// Update info panel of this project
	}
	
	
	// Adjust size of image to dimension of Display canvas.
	private synchronized boolean resizeImage(PvDisplayIf disp) 
	{
		if (m_currentlyResizing)
		{
			PsDebug.message("currently resizing!");
			return false;
		}
		else
			m_currentlyResizing = true;

		// PsDebug.message("resize");

		if (disp == null)
		{
			notifyAll();
			m_currentlyResizing = false;
			return false;
		}

		boolean resized	= false;

		if (disp == m_disp) 
		{
			Dimension dim = m_disp.getSize();
			if (dim.height>0 && dim.width>0 &&
				 (m_mis==null || m_imageHeight != dim.height || m_imageWidth != dim.width)) 
			{
				// Set resolution
				m_imageHeight = dim.height;
				m_imageWidth = dim.width;
				int oldNumBlocksX = m_numBlocksX;
				int oldNumBlocksY = m_numBlocksY;
				m_numBlocksX = (int) Math.ceil((double)(m_imageWidth) / m_blockSize.getValue());
				m_numBlocksY = (int) Math.ceil((double)(m_imageHeight) / m_blockSize.getValue());

				// Set up fluid solver
				if (m_fluidSolver.size == 0)
					initFluidSolver();
				else if (m_fluidSolver.n != m_numBlocksX - 2 || m_fluidSolver.m != m_numBlocksY - 2)
				{
					// PsDebug.message("resize fluidSolver");
					m_fluidSolver.resizeArray(m_numBlocksX - 2, m_numBlocksY - 2);

					// Clone oldFluidSolver
					try { m_oldFluidSolver = m_fluidSolver.clone(); }
					catch (CloneNotSupportedException e) { PsDebug.warning("Clone of fluidSolver not supported!"); }
				}
				
				// Resize mouse input to current size m_numBlocksX * m_numBlocksY in a convenient way
				resizeOldInputDensity(oldNumBlocksX, oldNumBlocksY);
				
				// Reset animation
				initAnimation();
				startAnimation();
				
				// Reset mouse traces
				resetMouseTraces();
				
				// Reset canvas, data field
				m_density = new PdVector(m_imageWidth*m_imageHeight);
				if(colorON) {
					m_density2 = new PdVector(m_imageWidth*m_imageHeight);
					m_density3 = new PdVector(m_imageWidth*m_imageHeight);
				}
				m_pix = new PiVector(m_imageWidth*m_imageHeight);
				// Create new MemoryImageSource
				m_mis = new MemoryImageSource(m_imageWidth, m_imageHeight, m_pix.m_data, 0, m_imageWidth);
				// Enable the possibility to update the pixels in m_pix
				m_mis.setAnimated(true);
				// Create image which will be painted as background in PvDisplay
				m_image = ((Component)m_disp).createImage(m_mis);
				m_disp.setBackgroundImage(m_image);
				m_disp.update(null);
				
				resized	= true;
			}

			notifyAll();
			m_currentlyResizing = false;
			return resized;
		}
		notifyAll();
		m_currentlyResizing = false;
		return false;
	}

	/** Invoked when component has been shown. */
	public void componentShown(ComponentEvent comp) {}
	/** Invoked when component has been hidden. */
	public void componentHidden(ComponentEvent comp) {}
	/** Invoked when component has been moved. */
	public void componentMoved(ComponentEvent comp) {}
	/** When component has been resized all images must be resized. */
	public void componentResized(ComponentEvent comp) 
	{
		// PsDebug.message("componentResized");
		// Adjust sizes of images to dimension of display canvas
		Object source = comp.getSource();
		if (source == m_disp) {
			if (!resizeImage(m_disp))
				return;
			m_disp.update(null);
		}
	}

	// Get display
	public PvDisplayIf getDisp() 
	{
		// PsDebug.message("getDisp");
		if (m_disp != null)
			return m_disp;
		
		if (getDisplay() != null) {
			// If available then use default project display.
			m_disp = getDisplay();
		} else {
			// Get viewer and ask for another, new display
			PvViewerIf viewer = getViewer();
			// Create right window and add a clone of the torus geometry.
			m_disp = viewer.newDisplay("Stable Fluids", false);
		}
		m_disp.setBackgroundColor(Color.white);
		m_disp.addPickListener(this);
//		((PvCameraIf)m_disp).addCameraListener(m_disp.getCamera());
		// PsDebug.message("adding component listener...");
		((Component)m_disp).addComponentListener(this);
		// PsDebug.message("added component listener!");
		m_disp.addCameraListener(new PvCameraListenerIf()
		{
			@Override
			public void pickCamera(PvCameraEvent cameraEvent) {/*PsDebug.message("pickCamera");*/}
			@Override
			public String getName() {return null;}
			@Override
			public void dragCamera(PvCameraEvent cameraEvent){painterDragCamera(cameraEvent);}
		});

		return m_disp;
	}

	// Compute image pixel values
	private synchronized void computeImage()
	{
		if(colorChange) { colorChange = false;} else {
			// PsDebug.message("Compute image!");
			// PsDebug.initTime();
	
			// Add forces on corresponding mouse trace
			if (m_forceTraceX.getSize() == m_forceTraceY.getSize())
				for (int i=0; i<m_forceTraceX.getSize()-1; i++)
					addForce(m_forceTraceX.getEntry(i), m_forceTraceY.getEntry(i), m_forceTraceX.getEntry(i+1), m_forceTraceY.getEntry(i+1), m_forceRadius.getValue());
	
			// Add densities on corresponding mouse trace
			addDensityTrace();
			
			// Delete mouse traces
			resetMouseTraces();
			
			// Solve fluid
			if(isFrozen) {} else {
				m_fluidSolver.velocitySolver();
				m_fluidSolver.densitySolver();
			}
	
			try { m_oldFluidSolver = m_fluidSolver.clone(); }
			catch (CloneNotSupportedException e) { PsDebug.warning("Clone of fluidSolver not supported!"); }
			
	        // Copy data from fluid solver into own data array... 
			boolean fluidSolverChanged = false;
			int fsWidth, fsHeight, blockSize, blockRemain;
			for (int y=0; y<m_imageHeight; y+=m_blockSize.getValue())
			{
				for (int x=0; x<m_imageWidth; x+=m_blockSize.getValue())
				{
					blockSize = m_blockSize.getValue();
					fsWidth = m_fluidSolver.n + 2;
					fsHeight = m_fluidSolver.m + 2;
	
					// If blockSize or canvas size was changed in the user interface, we first need to wait,
					// until size of fluidSolver is adjusted in changeBlockSize or resizeImage
					while (! (blockSize == m_oldBlockSize.getValue()
							&& fsWidth * blockSize >= m_imageWidth
							&& fsHeight * blockSize >= m_imageHeight
							&& (fsWidth - 1) * blockSize <= m_imageWidth
							&& (fsHeight - 1) * blockSize <= m_imageHeight)
							&& m_density.getSize() == m_imageWidth*m_imageHeight
							&& m_density2.getSize() == m_imageWidth*m_imageHeight
							&& m_density3.getSize() == m_imageWidth*m_imageHeight)
					{
						try { wait(); }
						catch (InterruptedException e) { PsDebug.warning("InterruptedExpeption in computeImage!"); }
						fluidSolverChanged = true;
						blockSize = m_blockSize.getValue();
						fsWidth = m_fluidSolver.n + 2;
						fsHeight = m_fluidSolver.m + 2;
					}
					// After changing the size of the fluidSolver, start double-loop from beginning
					if (fluidSolverChanged)
					{
						fluidSolverChanged = false;
						x = -m_blockSize.getValue();
						y = 0;
						continue;
					}
					
					// Otherwise, we can copy contend from fluidSolver into our pixelwise density-Array
					blockRemain = Math.min(m_blockSize.getValue(), m_imageWidth-x);
					for (int k=0; k<blockRemain; k++)
					{
						m_density.setEntry(I(x+k,y), Math.max(0, m_fluidSolver.d[Id(block(x), block(y))]));
						if(colorON) {
							m_density2.setEntry(I(x+k,y), Math.max(0, m_fluidSolver.d2[Id(block(x), block(y))]));
							m_density3.setEntry(I(x+k,y), Math.max(0, m_fluidSolver.d3[Id(block(x), block(y))]));
						}
					}
				}
	
				blockSize = m_blockSize.getValue();
				fsWidth = m_fluidSolver.n + 2;
				fsHeight = m_fluidSolver.m + 2;
	
				// If blockSize or canvas size was changed in the user interface, we first need to wait,
				// until size of fluidSolver is adjusted in changeBlockSize or resizeImage
				while (! (blockSize == m_oldBlockSize.getValue()
						&& fsWidth * blockSize >= m_imageWidth
						&& fsHeight * blockSize >= m_imageHeight
						&& (fsWidth - 1) * blockSize <= m_imageWidth
						&& (fsHeight - 1) * blockSize <= m_imageHeight)
						&& m_density.getSize() == m_imageWidth*m_imageHeight
						&& m_density2.getSize() == m_imageWidth*m_imageHeight
						&& m_density3.getSize() == m_imageWidth*m_imageHeight)
				{
					try { wait(); }
					catch (InterruptedException e) { PsDebug.warning("InterruptedExpeption in computeImage!"); }
					fluidSolverChanged = true;
					blockSize = m_blockSize.getValue();
					fsWidth = m_fluidSolver.n + 2;
					fsHeight = m_fluidSolver.m + 2;
				}
				if (fluidSolverChanged)
				{
					fluidSolverChanged = false;
					y = -m_blockSize.getValue();
					continue;
				}
				// Duplicate row column to fill blocks - only effective if blockSize > 1
				blockRemain = Math.min(m_blockSize.getValue(), m_imageHeight-y)-1;
				for (int k=1; k<=blockRemain; k++)
				{
					try {System.arraycopy(m_density.m_data, I(0,y), m_density.m_data, I(0,y+k), m_imageWidth);}
					catch (ArrayIndexOutOfBoundsException e) {}
					
					fluidSolverChanged = false;

					// Duplicate row column to fill blocks - only effective if blockSize > 1
					while (! (m_density.getSize() == m_imageWidth*m_imageHeight))
					{
						try { wait(); }
						catch (InterruptedException e) { PsDebug.warning("InterruptedExpeption in computeImage!"); }
						fluidSolverChanged = true;
					}
					if (fluidSolverChanged)
					{
						fluidSolverChanged = false;
						y = -m_blockSize.getValue();
						continue;
					}
					blockRemain = Math.min(m_blockSize.getValue(), m_imageHeight-y)-1;
					for (k=1; k<=blockRemain; k++)
					{
						try {System.arraycopy(m_density.m_data, I(0,y), m_density.m_data, I(0,y+k), m_imageWidth);}
						catch (ArrayIndexOutOfBoundsException e) {}
					}
					
					if(colorON) {
						fluidSolverChanged = false;
						// Duplicate row column to fill blocks - only effective if blockSize > 1
						while (! (m_density2.getSize() == m_imageWidth*m_imageHeight))
						{
							try { wait(); }
							catch (InterruptedException e) { PsDebug.warning("InterruptedExpeption in computeImage!"); }
							fluidSolverChanged = true;
						}
						if (fluidSolverChanged)
						{
							fluidSolverChanged = false;
							y = -m_blockSize.getValue();
							continue;
						}
						blockRemain = Math.min(m_blockSize.getValue(), m_imageHeight-y)-1;
						for (k=1; k<=blockRemain; k++)
						{
							try {System.arraycopy(m_density2.m_data, I(0,y), m_density2.m_data, I(0,y+k), m_imageWidth);}
							catch (ArrayIndexOutOfBoundsException e) {}
						}
						
						fluidSolverChanged = false;

						// Duplicate row column to fill blocks - only effective if blockSize > 1
						while (! (m_density3.getSize() == m_imageWidth*m_imageHeight))
						{
							try { wait(); }
							catch (InterruptedException e) { PsDebug.warning("InterruptedExpeption in computeImage!"); }
							fluidSolverChanged = true;
						}
						if (fluidSolverChanged)
						{
							fluidSolverChanged = false;
							y = -m_blockSize.getValue();
							continue;
						}
						blockRemain = Math.min(m_blockSize.getValue(), m_imageHeight-y)-1;
						for (k=1; k<=blockRemain; k++)
						{
							try {System.arraycopy(m_density3.m_data, I(0,y), m_density3.m_data, I(0,y+k), m_imageWidth);}
							catch (ArrayIndexOutOfBoundsException e) {}
						}
					}
				}
			}
			// ...and finally canvas
			computeColors();
			m_mis.newPixels(0, 0, m_imageWidth, m_imageHeight);
	
			// PsDebug.message("Seconds per Frame: " + PsDebug.getTimeUsed());
		}
	}
	

	// Compute color array from an array of scalar integer values
	private void computeColors()
	{	
		int red, green, blue;
		for (int x=0; x<m_imageWidth; x++)
		{
			for (int y=0; y<m_imageHeight; y++)
			{
				if(colorON) {
					if ((int)Math.round(m_density.getEntry(I(x,y))*255) > 255) {
						red = 0;
					} else if ((int)Math.round(m_density.getEntry(I(x,y))*255) < 0) {
						red = 255;
					} else {
						red = (int)Math.round((1.0-m_density.getEntry(I(x,y)))*255);
					}
					
					if ((int)Math.round(m_density2.getEntry(I(x,y))*255) > 255) {
						green = 0;
					} else if ((int)Math.round(m_density2.getEntry(I(x,y))*255) < 0) {
						green = 255;
					} else {
						green = (int)Math.round((1.0-m_density2.getEntry(I(x,y)))*255);
					}
					
					if ((int)Math.round(m_density3.getEntry(I(x,y))*255) > 255) {
						blue = 0;
					} else if ((int)Math.round(m_density3.getEntry(I(x,y))*255) < 0) {
						blue = 255;
					} else {
						blue = (int)Math.round((1.0-m_density3.getEntry(I(x,y)))*255);
					}
					
					int[] check = new int[3];
					/**
					if(red < 50) {check[0] = (255-red2[Id(x,y)]);} else {check[0] = (red + (255-red2[Id(x,y)])/2 );}
					if(red < 50) {check[1] = (255-red2[Id(x,y)]);} else {check[0] = (red + (255-red2[Id(x,y)])/2 );}
					if(red < 50) {check[2] = (255-red2[Id(x,y)]);} else {check[0] = (red + (255-red2[Id(x,y)])/2 );}
					*/
					
					check[0] = red;
					check[1] = green;
					check[2] = blue;
					m_pix.setEntry(I(x,y), PdColor.getDimmedColor( PdColor.getColor(255, check), (double) 1.0f)  );
					//if(x == 1 && y ==1)
					//PsDebug.message("Rot:" + String.valueOf(check[0]) + "Grün" + String.valueOf(check[1]) + "Blau" + String.valueOf(check[2]) );
					//m_pix.setEntry(I(x,y), PdColor.getColor(255, 0, 0, 255) );
				} else {
					if ((int)Math.round(m_density.getEntry(I(x,y))*255) > 255) {
						m_pix.setEntry(I(x,y), PdColor.hsv2rgbAsInt(0, 0, 0));
					} else if ((int)Math.round(m_density.getEntry(I(x,y))*255) < 0) {
						m_pix.setEntry(I(x,y), PdColor.hsv2rgbAsInt(0, 0, 255));
					} else {
						m_pix.setEntry(I(x,y), PdColor.hsv2rgbAsInt(0, 0, (int)Math.round((1.0-m_density.getEntry(I(x,y)))*255)));
					}
				}
			}
		}
	}
	
	/**
	 * Method is called from display when a user picks into the display
	 * in initial-pick mode.
	 * @param		pickEvent		Pick event issued by the display on left mouse pick.
	 * @see			jv.project.PvPickListenerIf
	 */
	public void pickInitial(PvPickEvent pickEvent)
	{
		// PsDebug.message("pickInitial");
		if (pickEvent.getSource() == m_disp)
		{
			Point loc = pickEvent.getLocation();
			m_densityTraceX.addEntry(loc.x);
			m_densityTraceY.addEntry(loc.y);
			loc.x = Math.floorDiv(loc.x, m_blockSize.getValue()) * m_blockSize.getValue();
			loc.y = Math.floorDiv(loc.y, m_blockSize.getValue()) * m_blockSize.getValue();
		}
	}

	/**
	 * Method is called from display when a drags picks into the display in scale mode.
	 * @param		cameraEvent		Camera event issued by the display on right mouse drag.
	 * @see			jv.project.PvPickListenerIf
	 */
	private void painterDragCamera(PvCameraEvent cameraEvent)
	{
		// PsDebug.message("painterDragCamera");
		if (cameraEvent.getSource() == m_disp && m_disp.getMajorMode() == PvDisplayIf.MODE_SCALE)
		{
			Point loc = cameraEvent.getLocation();
			m_forceTraceX.addEntry(loc.x);
			m_forceTraceY.addEntry(loc.y);
			loc.x = Math.floorDiv(loc.x, m_blockSize.getValue()) * m_blockSize.getValue();
			loc.y = Math.floorDiv(loc.y, m_blockSize.getValue()) * m_blockSize.getValue();
		}
	}

	private synchronized boolean changeBlockSize()
	{
		// PsDebug.message("changeBlockSize");
		// PsDebug.message("change blockSize from " + String.valueOf(m_oldBlockSize.getValue()) + " to " + String.valueOf(m_blockSize.getValue()));

		// Change resolutions
		int oldNumBlocksX = m_numBlocksX;
		m_numBlocksX = (int) Math.ceil((double)(m_imageWidth) / m_blockSize.getValue());
		m_numBlocksY = (int) Math.ceil((double)(m_imageHeight) / m_blockSize.getValue());

		// Save current fluidSolver
		try { m_oldFluidSolver = m_fluidSolver.clone(); }
		catch (CloneNotSupportedException e) { PsDebug.warning("Clone of fluidSolver not supported!"); }

		// Change resolution from fluidSolver
		m_fluidSolver.setup(m_numBlocksX - 2, m_numBlocksY - 2, (float)m_dt.getValue(), colorON);
		
		// Save current oldInputDensity
    	PdVector oidOld = (PdVector) m_oldInputDensity.clone();

    	// Resize oldInputDensity-Array to current size
    	m_oldInputDensity = new PdVector(m_numBlocksX * m_numBlocksY);

		// Average colors from blocks
		int oldSize = m_oldBlockSize.getValue();
		int newSize = m_blockSize.getValue();
		double temp, addThis;
		for (int x = 0; x < m_imageWidth; ++x)
		{
			for (int y = 0; y < m_imageHeight; ++y)
			{
				// FluidSolver
				m_fluidSolver.d[Id(block(x), block(y))] += m_oldFluidSolver.d[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.u[Id(block(x), block(y))] += m_oldFluidSolver.u[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.v[Id(block(x), block(y))] += m_oldFluidSolver.v[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.dOld[Id(block(x), block(y))] += m_oldFluidSolver.dOld[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.uOld[Id(block(x), block(y))] += m_oldFluidSolver.uOld[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.vOld[Id(block(x), block(y))] += m_oldFluidSolver.vOld[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				if(colorON) {
					m_fluidSolver.d2[Id(block(x), block(y))] += m_oldFluidSolver.d2[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
					m_fluidSolver.d3[Id(block(x), block(y))] += m_oldFluidSolver.d3[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
					m_fluidSolver.d2Old[Id(block(x), block(y))] += m_oldFluidSolver.d2Old[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
					m_fluidSolver.d3Old[Id(block(x), block(y))] += m_oldFluidSolver.d3Old[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				}
				
				// Old input density
				temp = m_oldInputDensity.getEntry(Id(block(x), block(y)));
				addThis = oidOld.getEntry(Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX));
				m_oldInputDensity.setEntry(Id(block(x), block(y)), temp + addThis);
			}
		}
		for (int x = 0; x < m_numBlocksX; ++x)
		{
			for (int y = 0; y < m_numBlocksY; ++y)
			{
				// FluidSolver
				m_fluidSolver.d[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.u[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.v[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.dOld[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.uOld[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.vOld[Id(x, y)] /= (newSize*newSize);
				if(colorON) {
					m_fluidSolver.d2Old[Id(x, y)] /= (newSize*newSize);
					m_fluidSolver.d3Old[Id(x, y)] /= (newSize*newSize);
					m_fluidSolver.d2[Id(x, y)] /= (newSize*newSize);
					m_fluidSolver.d3[Id(x, y)] /= (newSize*newSize);
				}
				
				// Old input density
				m_oldInputDensity.setEntry(Id(x, y), m_oldInputDensity.getEntry(Id(x, y)) / (newSize * newSize));
			}
		}
		
		// Save current state
		m_oldBlockSize.setValue(m_blockSize.getValue());
		try { m_oldFluidSolver = m_fluidSolver.clone(); }
		catch (CloneNotSupportedException e) { PsDebug.warning("Clone of fluidSolver not supported!"); }
		
		// Notify currently waiting threads (e.g. in computeImage), that fluidSolver is ready to use
		notifyAll();
		
		return true;
	}

	
	// Init and configure fluidSolver for current imageWidth, imageHeight, default constants
	private void initFluidSolver()
	{
		m_fluidSolver = new FluidSolver();
		m_fluidSolver.setup(m_numBlocksX - 2, m_numBlocksY - 2, (float)m_dt.getValue(), colorON);

		// configure fluidSolver
		m_fluidSolver.setDiff((float) m_diffusion.getValue() / 100000);
		m_fluidSolver.setBuoyancy(10.0f * (float) m_buoyancy.getValue());
		m_fluidSolver.setVorticity((float) m_vorticity.getValue());
		m_fluidSolver.setVisc((float) m_viscosity.getValue());
		
		// Clone oldFluidSolver
		try { m_oldFluidSolver = m_fluidSolver.clone(); }
		catch (CloneNotSupportedException e) { PsDebug.warning("Clone of fluidSolver not supported!"); }
		
	}
	
    /**
     * Calculate the mouse input force for each cell
     * in the fluid grid. We add force linearly decreasing
     * in circles around the "old" mouse position (oldX, oldY).
     * If m_blockSize > 1, then in one block, the average is taken 
     * over all pixels in that block.
     * Furthermore, since there are possibly multiple forces added in one
     * timestep, we just add the maximum of all input forces over all mouse input
     * force points to every block.
     *
     * @param oldX: first pixel coordinate in canvas, around which force is added.
     * @param oldY: second pixel coordinate in canvas, around which force is added.
     * @param newX: force in x-direction is determined by pixel difference newX-oldX.
     * @param newY: force in y-direction is determined by pixel difference newY-oldY.
     * @param radius: force is added up to radius pixels around old mouse position.
     **/
    private void addForce(int oldX, int oldY, int newX, int newY, int radius)
    {
//    	boolean wasZero = true;
//    	PsDebug.message("addForce to x: " + String.valueOf(m_mouseX_Old) + ", y: " + String.valueOf(m_mouseY_Old));
//    	PsDebug.message("new x: " + String.valueOf(m_mouseX) + ", new y: " + String.valueOf(m_mouseY));

    	float[] lowerUbound = new float[m_numBlocksX*m_numBlocksY];
    	float[] lowerVbound = new float[m_numBlocksX*m_numBlocksY];
    	PdVector mousePos = new PdVector((double)oldX, (double)oldY);
    	
    	// Add forces to lower bound, linearly decreasing from center
		for (int x=Math.max(0, oldX-radius); x<Math.min(m_imageWidth, oldX+radius); ++x)
		{
			for (int y=Math.max(0, oldY-radius); y<Math.min(m_imageHeight, oldY+radius); ++y)
			{
//				if (m_fluidSolver.vOld[I(x, y)] != 0 || m_fluidSolver.uOld[I(x, y)] != 0)
//				{
//					if (wasZero)
//				    	PsDebug.message("at x: " + String.valueOf(x) + ", y: " + String.valueOf(y) + " vOld/uOld is not zero");
//					wasZero = false;
//				}
//		    	PsDebug.message("dist == " + String.valueOf(mousePos.dist(new PdVector((double)x, (double)y))));
				double lambda = Math.max(0.0, (radius - mousePos.dist(new PdVector((double)x, (double)y))) / (double)radius);
				
//				PsDebug.message("x: " + String.valueOf(x) + ", y: " + String.valueOf(y)
//							+ ", lambda: " + String.valueOf(lambda));
				
				lowerUbound[Id(block(x), block(y))] += (float)(m_forceConst.getValue() * lambda * (newX - oldX) / m_blockSize.getValue());
				lowerVbound[Id(block(x), block(y))] += (float)(m_forceConst.getValue() * lambda * (newY - oldY) / m_blockSize.getValue());
			}
		}

		// DIESE GRENZEN TESTEN
		int uSign = newX - oldX >= 0 ? 1 : -1, vSign = newY - oldY >= 0 ? 1 : -1;
		for (int x=Math.max(0, block(oldX-radius)); x<Math.min(m_numBlocksX, block(oldX+radius)); ++x)
		{
			for (int y=Math.max(0, block(oldY-radius)); y<Math.min(m_numBlocksY, block(oldY+radius)); ++y)
			{
				// Average
				lowerUbound[Id(x, y)] /= (m_blockSize.getValue() * m_blockSize.getValue());
				lowerVbound[Id(x, y)] /= (m_blockSize.getValue() * m_blockSize.getValue());
				
				// Put maximum of input force that is already there and new forcecone into fluidSolver
				m_fluidSolver.uOld[Id(x, y)] = uSign * Math.max(uSign * lowerUbound[Id(x, y)], uSign * m_fluidSolver.uOld[Id(x, y)]);
				m_fluidSolver.vOld[Id(x, y)] = vSign * Math.max(vSign * lowerVbound[Id(x, y)], vSign * m_fluidSolver.vOld[Id(x, y)]);
			}
		}
				
//		if (!wasZero)
//			PsDebug.warning("vOld/uOld war nicht null!");
    }
    
    /**
     * Add density splines according to the current mouse trace
     */
    private synchronized void addDensityTrace()
    {
    	// Security check
    	if (m_densityTraceX.getSize() <= 1 && m_densityTraceX.getSize() <= 1)
    	{
			PsDebug.warning("m_densityTraces do not have enough entries! return...");
			return;
    	}

    	PgBezierCurve bezier;
    	PdVector p0, p1, a, b, c, d;
    	PdVector pointAtT;

    	// Security check
		if (m_densityTraceX.getSize() != m_densityTraceY.getSize())
		{
    		PsDebug.warning("Density traces do not have same size! return...");
    		return;
		}

    	// For all: mouse points in the trace
		for (int k = 1; k <= m_densityTraceX.getSize() - 2; ++k)
		{
			// Step 1: Calculate some points
			a = new PdVector(m_densityTraceX.getEntry(k), m_densityTraceY.getEntry(k));
			b = new PdVector(m_densityTraceX.getEntry(k+1), m_densityTraceY.getEntry(k+1));
			c = new PdVector(m_densityTraceX.getEntry(k-1), m_densityTraceY.getEntry(k-1));
			d = new PdVector(2);
			if (k < m_densityTraceX.getSize() - 2)
				d = new PdVector(m_densityTraceX.getEntry(k+2), m_densityTraceY.getEntry(k+2));
			p0 = (PdVector) a.clone();
			p0.sub(c);
			p0.multScalar((double)(2.0)/5 * a.dist(b) / a.dist(c));
			p0.add(a);
			p1 = (PdVector) b.clone();
			p1.sub(d);
			p1.multScalar((double)(2.0)/5 * a.dist(b) / b.dist(d));
			p1.add(b);
			bezier = new PgBezierCurve(0);
			bezier.setDimOfVertices(2);

			// Step 2: Construct the bezier curve from the control points
			// If there is just one single new point without other information...
			if (a.getEntry(0) == -1 && a.getEntry(1) == -1)
			{
				// ...make control polygon with 1 control point
				bezier.setNumControlPoints(1);
				bezier.setControlPoint(0, b);
			}
			// Otherwise, if there is information where the curve comes from...
			else if (c.getEntry(0) != -1 && c.getEntry(1) != -1)
			{
				// ...and if there is information where the curve goes to...
				if (k < m_densityTraceX.getSize() - 2)
				{
					// ...make control polygon with 4 control points
					bezier.setNumControlPoints(4);
					bezier.setControlPoint(0, a);
					bezier.setControlPoint(1, p0);
					bezier.setControlPoint(2, p1);
					bezier.setControlPoint(3, b);
				}
				// ...but no information, where the curve goes to...
				else
				{
					// ...make control polygon with 3 control points
					bezier.setNumControlPoints(3);
					bezier.setControlPoint(0, a);
					bezier.setControlPoint(1, p0);
					bezier.setControlPoint(2, b);
				}
			}
			// ...but if there is no information, where the curve comes from...
			else
			{
				// ...and if there is information where the curve goes to...
				if (k < m_densityTraceX.getSize() - 2)
				{
					// ...make control polygon with 3 control points
					bezier.setNumControlPoints(3);
					bezier.setControlPoint(0, a);
					bezier.setControlPoint(1, p1);
					bezier.setControlPoint(2, b);
				}
				// ...and no information, where the curve goes to...
				else
				{
					// ...make control polygon with 2 control points
					bezier.setNumControlPoints(2);
					bezier.setControlPoint(0, a);
					bezier.setControlPoint(1, b);
				}
			}
			
			// Step 3: Add density along bezier curve
			pointAtT = new PdVector(2);
			for (int i = 0; i <= NUM_PARTS_DENSITY_TRACE_SPLINES; ++i)
			{
				bezier.eval(pointAtT, (double)(i)/NUM_PARTS_DENSITY_TRACE_SPLINES);
				addDensityPoint((int)Math.round(pointAtT.getEntry(0)), (int)Math.round(pointAtT.getEntry(1)), m_densityRadius.getValue());
			}
		}
		
		// Save this state of the mouse input density
		for (int x=0; x<m_numBlocksX; ++x)
		{
			for (int y=0; y<m_numBlocksY; ++y)
			{
				if (whichColor == 1 || whichColor == 5)
					m_oldInputDensity.setEntry(Id(x, y), m_fluidSolver.dOld[Id(x, y)]);
				else if (whichColor == 2)
					m_oldInputDensity.setEntry(Id(x, y), m_fluidSolver.d2Old[Id(x, y)]);
				else if (whichColor == 3 || whichColor == 4 || whichColor == 6 || whichColor == 7 || whichColor == 8)
					m_oldInputDensity.setEntry(Id(x, y), m_fluidSolver.d3Old[Id(x, y)]);
			}
		}
    }
    
    /**
     * Add density linearly decreasing in circles around the pixel position (pX, pY).
     *
     * @param pX: density is added in the block with first pixel coordinate pX.
     * @param pY: density is added in the block with second pixel coordinate pY.
     * @param radius: density is added up to radius around mouse position.
     **/
    private void addDensityPoint(int pX, int pY, int radius)
    {
//    	PsDebug.message("addForce to x: " + String.valueOf(m_mouseX_Old) + ", y: " + String.valueOf(m_mouseY_Old));

    	float[] lowerDbound = new float[m_numBlocksX*m_numBlocksY];
    	PdVector mousePos = new PdVector((double)pX, (double)pY);
    	final float densityNormConst = 3/(float)(radius*radius)/(float)Math.PI;
    	double temp;
    	
    	// Add densities to lower bound, linearly decreasing from center
		for (int x=Math.max(0, pX-radius); x<Math.min(m_imageWidth, pX+radius); x++)
		{
			for (int y=Math.max(0, pY-radius); y<Math.min(m_imageHeight, pY+radius); y++)
			{
				double lambda = Math.max(0.0, (radius - mousePos.dist(new PdVector((double)x, (double)y))) / (double)radius);
				lowerDbound[Id(block(x), block(y))] += (float)(lambda * m_densityConst.getValue() * densityNormConst);
			}
		}
		// Loop over all blocks in the square of the added density cone
		for (int x=Math.max(0, block(pX-radius)); x<Math.min(m_numBlocksX, block(pX+radius)); ++x)
		{
			for (int y=Math.max(0, block(pY-radius)); y<Math.min(m_numBlocksY, block(pY+radius)); ++y)
			{
				// Average
				lowerDbound[Id(x, y)] /= (m_blockSize.getValue() * m_blockSize.getValue());

				// Put maximum of input force that is already there, that was there last frame
				// and new density cone into fluidSolver
				if(whichColor == 1) {
					// m_fluidSolver.dOld[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.dOld[Id(x, y)]);
					temp = m_fluidSolver.dOld[Id(x,y)];
					m_fluidSolver.dOld[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
				} else if(whichColor == 2) {
					// m_fluidSolver.d2Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d2Old[Id(x, y)]);
					temp = m_fluidSolver.d2Old[Id(x,y)];
					m_fluidSolver.d2Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
				} else if(whichColor == 3) {
					// m_fluidSolver.d3Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d3Old[Id(x, y)]);
					temp = m_fluidSolver.d3Old[Id(x,y)];
					m_fluidSolver.d3Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
				} else if(whichColor == 4) {
					// m_fluidSolver.d3Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d3Old[Id(x, y)]);
					temp = m_fluidSolver.d3Old[Id(x,y)];
					m_fluidSolver.d3Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
					// m_fluidSolver.d2Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d2Old[Id(x, y)]);
					temp = m_fluidSolver.d2Old[Id(x,y)];
					m_fluidSolver.d2Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
				} else if(whichColor == 5) {
					// m_fluidSolver.dOld[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.dOld[Id(x, y)]);
					temp = m_fluidSolver.dOld[Id(x,y)];
					m_fluidSolver.dOld[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
					// m_fluidSolver.d2Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d2Old[Id(x, y)]);
					temp = m_fluidSolver.d2Old[Id(x,y)];
					m_fluidSolver.d2Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
				} else if(whichColor == 6) {
					// m_fluidSolver.dOld[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.dOld[Id(x, y)]);
					temp = m_fluidSolver.dOld[Id(x,y)];
					m_fluidSolver.dOld[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
					// m_fluidSolver.d3Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d3Old[Id(x, y)]);
					temp = m_fluidSolver.d3Old[Id(x,y)];
					m_fluidSolver.d3Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
				}
				// else if(whichColor == 7) {
				// 	// m_fluidSolver.d3Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d3Old[Id(x, y)]);
				// 	temp = m_fluidSolver.d3Old[Id(x,y)];
				// 	m_fluidSolver.d3Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
				// 	// m_fluidSolver.d2Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d2Old[Id(x, y)]);
				// 	temp = m_fluidSolver.d2Old[Id(x,y)];
				// 	m_fluidSolver.d2Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
				// 	// m_fluidSolver.dOld[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.dOld[Id(x, y)])/2;
				// 	temp = m_fluidSolver.dOld[Id(x,y)];
				// 	m_fluidSolver.dOld[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) (m_oldInputDensity.getEntry(Id(x, y))/2)))/2;
				// } else if(whichColor == 8) {
				// 	// m_fluidSolver.d3Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d3Old[Id(x, y)]);
				// 	temp = m_fluidSolver.d3Old[Id(x,y)];
				// 	m_fluidSolver.d3Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y)))); 
				// 	// m_fluidSolver.d2Old[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.d2Old[Id(x, y)])/2;
				// 	temp = m_fluidSolver.d2Old[Id(x,y)];
				// 	m_fluidSolver.d2Old[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y))/2))/2;
				// 	// m_fluidSolver.dOld[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.dOld[Id(x, y)])/2;
				// 	temp = m_fluidSolver.dOld[Id(x,y)];
				// 	m_fluidSolver.dOld[Id(x, y)] = (float) Math.max(temp, Math.max(0.0, Math.max(lowerDbound[Id(x, y)], temp) - (float) m_oldInputDensity.getEntry(Id(x, y))/2))/2;
				// }
			}
		}
    }
    
    /**
     * Resets the force- and densityTraces, i.e. deletes every entry except for 
     * the last one, which is then the next "old" mouse position.
     **/
    private void resetMouseTraces()
    {
    	// Force traces
		m_forceTraceX = new PiVector();
		m_forceTraceX.setSize(0);
		m_forceTraceY = new PiVector();
		m_forceTraceY.setSize(0);
		
		// Density traces - the last two entries are saved for being able to draw C2-splines.
		// If there is no new mouse point, set everything to (-1, -1).
		if (m_densityTraceX.getSize() == m_densityTraceY.getSize() && m_densityTraceX.getSize() > 2)
		{
			m_densityTraceX.setEntry(0, m_densityTraceX.getEntry(m_densityTraceX.getSize()-2));
			m_densityTraceX.setEntry(1, m_densityTraceX.getLastEntry());
			m_densityTraceX.setSize(2);
			m_densityTraceY.setEntry(0, m_densityTraceY.getEntry(m_densityTraceY.getSize()-2));
			m_densityTraceY.setEntry(1, m_densityTraceY.getLastEntry());
			m_densityTraceY.setSize(2);
		}
		else
		{
			if (m_densityTraceX.getSize() != m_densityTraceY.getSize())
				PsDebug.warning("m_densityTraceX and ..Y have different size! Reset them to {(-1,-1), (-1,-1)}...");

			m_densityTraceX = new PiVector();
			m_densityTraceX.setSize(2);
			m_densityTraceX.setEntry(0, -1);
			m_densityTraceX.setEntry(1, -1);
			m_densityTraceY = new PiVector();
			m_densityTraceY.setSize(2);
			m_densityTraceY.setEntry(0, -1);
			m_densityTraceY.setEntry(1, -1);
		}
	}
    
    // Cuts off or expands by some zeros the array m_oldInputDensity. Shall be invoked in resizeImage.
    private void resizeOldInputDensity(int oldNumBlocksX, int oldNumBlocksY)
    {
    	// Clone array
    	PdVector temp = (PdVector) m_oldInputDensity.clone();
    	
    	// Resize oldInputDensity-Array to current size
    	m_oldInputDensity = new PdVector(m_numBlocksX * m_numBlocksY);
    	
    	// Copy old values back
    	for (int x=0; x<Math.min(oldNumBlocksX, m_numBlocksX); ++x)
    		for (int y=0; y<Math.min(oldNumBlocksY, m_numBlocksY); ++y)
    			m_oldInputDensity.setEntry(Id(x,y), temp.getEntry(Id(x,y,oldNumBlocksX)));
    }
	

	/* Utility indexing functions */

    // Calculates the index of a 1d-Array for given x,y-coordinates of 2d *image*-array
	private int I(int x, int y)	{return x + (m_imageWidth*y);}
    // Calculates the index of a 1d-Array for given x,y-coordinates of 2d *datagrid*-array
	private int Id(int x, int y) {return x + (m_numBlocksX*y);}
	// Same as above, but for "custom" (parameter) width of data array
	private int Id(int x, int y, int arrayWidth) {return x + (arrayWidth*y);}
    // Calculates the index of a 1d-Array for given x,y-coordinates of 2d-array, which has a 1-pix-border around
	private int I2(int x, int y)	{return x + ((m_imageWidth+2)*y);}
	// "Inverse functions" of I
	private int ind2x(int ind) {return (ind % m_imageWidth);}
	private int ind2y(int ind) {return Math.floorDiv(ind, m_imageWidth);}
	// Calculates the block (starting from 0), for the current m_blockSize, that some given pixel pix is in
	private int block(int pix) {return Math.floorDiv(pix, m_blockSize.getValue());}
	// Same as above, but for "custom" (parameter) blockSize
	private int block(int pix, int blockSize) {return Math.floorDiv(pix, blockSize);}
	// Returns the corner pixel value for a given block index
	private int cornerPix(int block) {return block * m_blockSize.getValue();}
	// Same as above, but for "custom" (parameter) blockSize
	private int cornerPix(int block, int blockSize) {return block * blockSize;}
	
	public PuDouble getTime() { return m_time; }
}