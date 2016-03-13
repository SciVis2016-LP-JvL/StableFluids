package StableFluids;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
// import java.awt.GridLayout;
import java.awt.Image;
//import java.awt.Label;
// import java.awt.Panel;
import java.awt.Point;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.image.MemoryImageSource;
//import java.awt.event.MouseEvent;
//import java.awt.event.MouseListener;

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

import jvx.curve.PgBezierCurve;

/**
 * Project visualizing smoke like [Jos Stam 1999].
 * 
 * @author		Lukas Polthier, Johannes von Lindheim, Konrad Polthier (Julia Set project)
 * @version		0.0, 20.12.15 reused JuliaSets project as template for Stable Fluids application
 */
@SuppressWarnings("serial")
public class Painter extends PjProject implements ComponentListener {
	
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

	// Handling mouse input
	private		PiVector			m_forceTraceX;
	private		PiVector			m_forceTraceY;
	protected	PuInteger			m_forceRadius;
	protected	PuDouble			m_forceConst;
	private		PiVector			m_densityTraceX;
	private		PiVector			m_densityTraceY;
	protected	PuInteger			m_densityRadius;
	protected	PuDouble			m_densityConst;

	/**
	 * Determines size of uniformly colored pixels blocks in images.
	 * Using discr==1 leads to highest resolution where each image pixel is really computed.
	 */
	protected	PuInteger			m_blockSize;
	protected	PuInteger			m_oldBlockSize;
	
	private		PuDouble			m_time;
//	protected	Label				m_lTime;
	
	// Calculates the simulation data
	private 	FluidSolver			m_fluidSolver;
	private 	FluidSolver			m_oldFluidSolver;

    private 	PuDouble 			m_dt;
    
    private		PgBezierCurve		m_bezier;		
    
	public Painter() {
		super("Stable Fluids");
		myReset();
	}
	public void myReset()
	{
		// PsDebug.message("reset");
		m_pix					= new PiVector();
		m_density				= new PdVector();
		m_forceTraceX			= new PiVector();
		m_forceTraceY			= new PiVector();
		m_forceRadius			= new PuInteger("Force radius");
		m_forceConst			= new PuDouble("Force constant");
		m_densityTraceX			= new PiVector();
		m_densityTraceY			= new PiVector();
		m_densityRadius			= new PuInteger("Density radius");
		m_densityConst			= new PuDouble("Density constant");
		m_fluidSolver			= new FluidSolver();
		m_oldFluidSolver		= new FluidSolver();
		// // TEST
		// m_fluidSolver.setup(3, 3, (float)0.2);
		// m_fluidSolver.d[1] = (float) 0.5;
		// m_oldFluidSolver.setup(2, 3, (float)0.2);
		// PsDebug.message(String.valueOf(m_oldFluidSolver.d[1]));
		// m_oldFluidSolver.n = 1;
		// PsDebug.message(String.valueOf(m_oldFluidSolver.n));
		// try { m_oldFluidSolver = m_fluidSolver.clone(); }
		// catch (CloneNotSupportedException e) { PsDebug.warning("Clone of fluidSolver not supported!"); }
		// PsDebug.message(String.valueOf(m_oldFluidSolver.n));
		// PsDebug.message(String.valueOf(m_fluidSolver.n));
		// PsDebug.message(String.valueOf(m_oldFluidSolver.d[1]));
		// m_fluidSolver.d[1] = (float) 5.0;
		// PsDebug.message(String.valueOf(m_oldFluidSolver.d[1]));
		// m_oldFluidSolver.d[1] = (float) 3.0;
		// PsDebug.message(String.valueOf(m_oldFluidSolver.d[1]));
		// PsDebug.message(String.valueOf(m_fluidSolver.d[1]));
		// // TEST ENDE
		m_dt					= new PuDouble("Timestep");
		m_time					= new PuDouble("Time");
		m_blockSize				= new PuInteger("Block Size", this);
		m_oldBlockSize			= new PuInteger("Old block size", this);

		// For testing purposes
		m_bezier				= new PgBezierCurve(0);
		
		// Initialize animation
		initAnimation();

		if (getClass() == Painter.class) {
			init();
		}
	}

	// Called, when reset-Button is pressed
	public void buttonReset()
	{
		init();
		initAnimation();
		myStart();
	}

	public void init()
	{
		super.init();

		m_image					= null;
		m_pix.setSize(0);
		m_density.setSize(0);
		m_forceTraceX.setSize(0);
		m_forceTraceY.setSize(0);
		m_forceRadius.setBounds(1, 40, 1, 5);
		m_forceRadius.setValue(20);
		m_forceConst.setBounds(0.0, 0.6, 0.01, 0.1);
		m_forceConst.setValue(0.3);
		m_densityTraceX.setSize(0);
		m_densityTraceY.setSize(0);
		m_densityRadius.setBounds(1, 20, 1, 5);
		m_densityRadius.setValue(10);
		m_densityConst.setBounds(0.0, 600.0, 10.0, 100.0);
		m_densityConst.setValue(300.0);
		m_time.setValue(0);
		
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
		
		m_bezier.setDimOfVertices(2);
		m_bezier.setNumControlPoints(4);
		m_bezier.setControlPoint(0, new PdVector(5.0,5.0));
		m_bezier.setControlPoint(1, new PdVector(200,100));
		m_bezier.setControlPoint(2, new PdVector(180, 200));
		m_bezier.setControlPoint(3, new PdVector(50, 220));
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
		//m_disp.showScenegraph(false);
		
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
			// waitUntilNotWriting();
			changeBlockSize();

			computeImage();
			// m_disp.update(null);
			
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
		anim.setTime(0.0);
		anim.setTimeInterval(0., .2, .1, .1);

		// try 
		// {
			// Panel closeButtonPanel = (Panel) anim.getAnimationPanel().getComponent(1);
			// closeButtonPanel.remove(0);
			// anim.getAnimationPanel().remove(0);
			// Panel pArea = new Panel();
			// pArea.setLayout(new GridLayout(1, 2));
			// pArea.add(new Label("Time"));
			// m_lTime = new Label();
			// pArea.add(m_lTime);
			// anim.getAnimationPanel().add(pArea, new FlowLayout(FlowLayout.LEADING));
		// } catch(ArrayIndexOutOfBoundsException e){}
		
		// anim.stop();
		// anim.getAnimationPanel().setVisible(false);
		// anim.getAnimationPanel().setLocation(100, 350);
		// anim.setAnimationPanel(null);
	}
	
	// Start animation, if there is one
	private void startAnimation()
	{
		if (super.hasAnimation())
			super.getAnimation().start();
	}
	
	public boolean setTime(PsTimeEvent timeEvent) {
		// PsDebug.message("setTime");
		m_time.setValue(Math.round(timeEvent.getTime()));
		computeImage();
		m_disp.update(null);
//		m_lTime.setText(String.valueOf(m_time));
		if (super.hasAnimation())
			super.getAnimation().setTimeInterval(0., super.getAnimation().getMaxTime() + 1.);
		return super.update(this);				// Update info panel of this project
	}
	
	
	// Adjust size of image to dimension of Display canvas.
	private synchronized boolean resizeImage(PvDisplayIf disp) 
	{
		PsDebug.message("resize");

		if (disp == null)
		{
			notifyAll();
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
				m_numBlocksX = (int) Math.ceil((double)(m_imageWidth) / m_blockSize.getValue());
				m_numBlocksY = (int) Math.ceil((double)(m_imageHeight) / m_blockSize.getValue());

				// Set up fluid solver
				m_dt.setValue(0.2f);
				if (m_fluidSolver.size == 0)
				{
					PsDebug.message("create new fluidSolver");
					m_fluidSolver = new FluidSolver();
					m_fluidSolver.setup(m_numBlocksX - 2, m_numBlocksY - 2, (float)m_dt.getValue());
					m_oldFluidSolver = new FluidSolver();
					m_oldFluidSolver.setup(m_numBlocksX - 2, m_numBlocksY - 2, (float)m_dt.getValue());
				}
				else if (m_fluidSolver.n != m_numBlocksX || m_fluidSolver.m != m_numBlocksY)
				{
					// PsDebug.message("resize fluidSolver");
					// m_fluidSolver.resizeArray(m_numBlocksX - 2, m_numBlocksY - 2);
					// m_oldFluidSolver.resizeArray(m_numBlocksX - 2, m_numBlocksY - 2);
					PsDebug.message("create new fluidSolver");
					m_fluidSolver = new FluidSolver();
					m_fluidSolver.setup(m_numBlocksX - 2, m_numBlocksY - 2, (float)m_dt.getValue());
					// m_oldFluidSolver = new FluidSolver();
					// m_oldFluidSolver.setup(m_numBlocksX - 2, m_numBlocksY - 2, (float)m_dt.getValue());
				}
				
				// Reset animation
				initAnimation();
				startAnimation();
				
				// Reset mouse traces
				resetMouseTraces();
				
				// Reset canvas, data field
				m_density = new PdVector(m_imageWidth*m_imageHeight);
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
			return resized;
		}
		notifyAll();
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
			public void pickCamera(PvCameraEvent cameraEvent) {}
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
		// PsDebug.message("Compute image!");
		// PsDebug.initTime();

		// Add forces on corresponding mouse trace
		if (m_forceTraceX.getSize() == m_forceTraceY.getSize())
			for (int i=0; i<m_forceTraceX.getSize()-1; i++)
				addForce(m_forceTraceX.getEntry(i), m_forceTraceY.getEntry(i), m_forceTraceX.getEntry(i+1), m_forceTraceY.getEntry(i+1), m_forceRadius.getValue());

		// Add densities on corresponding mouse trace
		if (m_densityTraceX.getSize() == m_densityTraceY.getSize())
			for (int i=0; i<m_densityTraceX.getSize(); i++)
				addDensity(m_densityTraceX.getEntry(i), m_densityTraceY.getEntry(i), m_densityRadius.getValue());
		else
			PsDebug.warning("m_densityTraceX and ..Y have different size!");
		
		// Delete mouse traces
		resetMouseTraces();
		
		// Solve fluid
        m_fluidSolver.velocitySolver();
        m_fluidSolver.densitySolver();

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

				// If blockSize or canvas size were changed in the user interface,
				// we first need to wait, until size of fluidSolver is adjusted in changeBlockSize or resizeImage
				while (! (blockSize == m_oldBlockSize.getValue()
						&& fsWidth * blockSize >= m_imageWidth
						&& fsHeight * blockSize >= m_imageHeight
						&& (fsWidth - 1) * blockSize <= m_imageWidth
						&& (fsHeight - 1) * blockSize <= m_imageHeight)
						&& m_density.getSize() == m_imageWidth*m_imageHeight)
				{
					// PsDebug.message(String.valueOf(m_fluidSolver.size));
					// PsDebug.message(String.valueOf(fsWidth));
					// PsDebug.message(String.valueOf(fsHeight));
					// PsDebug.message(String.valueOf(m_imageWidth));
					// PsDebug.message(String.valueOf(m_imageHeight));
					// PsDebug.message(String.valueOf(m_blockSize.getValue()));
					// PsDebug.message(String.valueOf(m_oldBlockSize.getValue()));
					try { PsDebug.message("Wait"); wait(); }
					catch (InterruptedException e) { PsDebug.warning("InterruptedExpeption in computeImage!"); }
					fluidSolverChanged = true;
					blockSize = m_blockSize.getValue();
					fsWidth = m_fluidSolver.n + 2;
					fsHeight = m_fluidSolver.m + 2;
				}
				// PsDebug.message("done with waiting");
				// After changing the sized of the fluidSolver, start double-loop from beginning
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
					// try { m_density.setEntry(I(x+k,y), Math.max(0, m_fluidSolver.d[Id(block(x), block(y))])); }
					// catch(ArrayIndexOutOfBoundsException e)
					// {
					// 	PsDebug.message(String.valueOf(Id(block(x), block(y))) + " goes out of bounds!");
					// 	PsDebug.message("oldFluidSolverSize: " + String.valueOf(m_oldFluidSolver.size));
					// 	PsDebug.message("fluidSolverSize: " + String.valueOf(m_fluidSolver.size));
					// 	PsDebug.message("blockSize: " + String.valueOf(m_blockSize.getValue()));
					// 	PsDebug.message("oldBlockSize: " + String.valueOf(m_oldBlockSize.getValue()));
					// }
					m_density.setEntry(I(x+k,y), Math.max(0, m_fluidSolver.d[Id(block(x), block(y))]));
				}
			}
			// Duplicate row column to fill blocks - only effective if blockSize > 1
			while (! (m_density.getSize() == m_imageWidth*m_imageHeight))
			{
				PsDebug.message("in while");
				PsDebug.message(String.valueOf(m_density.getSize()));
				PsDebug.message(String.valueOf(m_fluidSolver.size));
				PsDebug.message(String.valueOf(m_imageWidth));
				PsDebug.message(String.valueOf(m_imageHeight));
				PsDebug.message(String.valueOf(m_blockSize.getValue()));
				PsDebug.message(String.valueOf(m_oldBlockSize.getValue()));
				try { PsDebug.message("Wait later"); wait(); }
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
			for (int k=1; k<=blockRemain; k++)
			{
				try {System.arraycopy(m_density.m_data, I(0,y), m_density.m_data, I(0,y+k), m_imageWidth);}
				catch (ArrayIndexOutOfBoundsException e)
				{
					PsDebug.message("I am out of bounds at arraycopy!");
					PsDebug.message(String.valueOf(m_density.getSize()));
					PsDebug.message(String.valueOf(I(0,y)));
					PsDebug.message(String.valueOf(I(0,y+k)));
					PsDebug.message(String.valueOf(m_fluidSolver.size));
					PsDebug.message(String.valueOf(m_imageWidth));
					PsDebug.message(String.valueOf(m_imageHeight));
					PsDebug.message(String.valueOf(m_blockSize.getValue()));
					PsDebug.message(String.valueOf(m_oldBlockSize.getValue()));
				}
			}
		}
		// ...and finally canvas
		computeColors();
		m_mis.newPixels(0, 0, m_imageWidth, m_imageHeight);

		// PsDebug.message("Seconds per Frame: " + PsDebug.getTimeUsed());
	}
	

	// Compute color array from an array of scalar integer values
	private void computeColors()
	{
		// Compute greyscale values according to density
		for (int x=0; x<m_imageWidth; x++)
			for (int y=0; y<m_imageHeight; y++)
				if ((int)Math.round(m_density.getEntry(I(x,y))*255) > 255)
//					PsDebug.warning("Color out of bounds!");
					m_pix.setEntry(I(x,y), PdColor.hsv2rgbAsInt(0, 0, 0));
				else if ((int)Math.round(m_density.getEntry(I(x,y))*255) < 0)
					m_pix.setEntry(I(x,y), PdColor.hsv2rgbAsInt(0, 0, 255));
				else
					m_pix.setEntry(I(x,y), PdColor.hsv2rgbAsInt(0, 0, (int)Math.round((1.0-m_density.getEntry(I(x,y)))*255)));
//					m_pix.setEntry(I(x,y), PdColor.hsv2rgbAsInt(0, 0, 1*255));
		
		// For testing purposes, show bezier curve
		final int noEvalsPerCtrlPt = 1000;
		PdVector pointAtT = new PdVector(2);
		pointAtT.setEntry(0, 15.0);
		pointAtT.setEntry(1, 15.0);
		for (int j=0; j<noEvalsPerCtrlPt; ++j)
		{
			// PsDebug.message("no control points == " + String.valueOf(m_bezier.getNumControlPoints()));
			// PsDebug.message("m_dim == " + String.valueOf(m_bezier.getDimOfVertices()));
			// PsDebug.message("Try to evaluate for t=" + String.valueOf((double)(j)/noEvalsPerCtrlPt) + "...");
			// for (int k=0; k<m_bezier.getNumControlPoints(); ++k)
				// PsDebug.message("x: " + String.valueOf((int)Math.round(m_bezier.getControlPoint()[k].getEntry(0))) + ", y: " + String.valueOf((int)Math.round(m_bezier.getControlPoint()[k].getEntry(1))));

			m_bezier.eval(pointAtT, (double)(j)/noEvalsPerCtrlPt);
			// PsDebug.message("x: " + String.valueOf((int)Math.round(pointAtT.getEntry(0))) + ", y: " + String.valueOf((int)Math.round(pointAtT.getEntry(1))));

			if (pointAtT.getSize() != 0)
				m_pix.setEntry(I((int)Math.round(pointAtT.getEntry(0)), (int)Math.round(pointAtT.getEntry(1))), PdColor.hsv2rgbAsInt(0, 1*255, 1*255));
			else
			{
				PsDebug.message("Size == 0. No. of control points: " + String.valueOf(m_bezier.getNumControlPoints()));
				for (int k=0; k<m_bezier.getNumControlPoints(); ++k)
					m_pix.setEntry(I((int)Math.round(m_bezier.getControlPoint()[k].getEntry(0)), (int)Math.round(m_bezier.getControlPoint()[k].getEntry(1))), PdColor.hsv2rgbAsInt(0, 1*255, 1*255));
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
		if (pickEvent.getSource() == m_disp)
		{
			Point loc = pickEvent.getLocation();
			m_densityTraceX.addEntry(loc.x);
			m_densityTraceY.addEntry(loc.y);
//			m_mouseX = loc.x;
//			m_mouseY = loc.y;
			loc.x = Math.floorDiv(loc.x, m_blockSize.getValue()) * m_blockSize.getValue();
			loc.y = Math.floorDiv(loc.y, m_blockSize.getValue()) * m_blockSize.getValue();
//	    	PsDebug.message("pickInitial x: " + String.valueOf(m_mouseX) + ", y: " + String.valueOf(m_mouseY));
			
//			PsDebug.message("dOld size:" + String.valueOf(m_fluidSolver.dOld.length));
//			m_fluidSolver.dOld[I(loc.x, loc.y)] = 100;
//			addDensity(10);
		}
	}
	/**
	 * Method is called from display when a drags picks into the display in scale mode.
	 * @param		cameraEvent		Camera event issued by the display on right mouse drag.
	 * @see			jv.project.PvPickListenerIf
	 */
	private void painterDragCamera(PvCameraEvent cameraEvent)
	{
		if (cameraEvent.getSource() == m_disp && m_disp.getMajorMode() == PvDisplayIf.MODE_SCALE)
		{
			Point loc = cameraEvent.getLocation();
			m_forceTraceX.addEntry(loc.x);
			m_forceTraceY.addEntry(loc.y);
//			m_mouseX = loc.x;
//			m_mouseY = loc.y;
			loc.x = Math.floorDiv(loc.x, m_blockSize.getValue()) * m_blockSize.getValue();
			loc.y = Math.floorDiv(loc.y, m_blockSize.getValue()) * m_blockSize.getValue();
//	    	PsDebug.message("dragCamera x: " + String.valueOf(m_mouseX) + ", y: " + String.valueOf(m_mouseY));
		}
	}

	private synchronized boolean changeBlockSize()
	{
		// PsDebug.message("changeBlockSize");
		PsDebug.message("change blockSize from " + String.valueOf(m_oldBlockSize.getValue()) + " to " + String.valueOf(m_blockSize.getValue()));

		// Change resolutions
		int oldNumBlocksX = m_numBlocksX;
		m_numBlocksX = (int) Math.ceil((double)(m_imageWidth) / m_blockSize.getValue());
		m_numBlocksY = (int) Math.ceil((double)(m_imageHeight) / m_blockSize.getValue());

		// BE CAREFUL NOT TO OVERWRITE IT WITH SOMETHING I ALREADY DID IN RESIZEIMAGE!
		// Save current fluidSolver
		try { m_oldFluidSolver = m_fluidSolver.clone(); }
		catch (CloneNotSupportedException e) { PsDebug.warning("Clone of fluidSolver not supported!"); }

		// Change resolution from fluidSolver
		m_fluidSolver.setup(m_numBlocksX - 2, m_numBlocksY - 2, (float)m_dt.getValue());

		// testoutputs
		// PsDebug.message("oldFluidSolverSize: " + String.valueOf(m_oldFluidSolver.size));
		// PsDebug.message("fluidSolverSize: " + String.valueOf(m_fluidSolver.size));

		// Average colors from blocks
		int oldSize = m_oldBlockSize.getValue();
		int newSize = m_blockSize.getValue();
		for (int x = 0; x < m_imageWidth; ++x)
		{
			for (int y = 0; y < m_imageHeight; ++y)
			{
				m_fluidSolver.d[Id(block(x), block(y))] += m_oldFluidSolver.d[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.u[Id(block(x), block(y))] += m_oldFluidSolver.u[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.v[Id(block(x), block(y))] += m_oldFluidSolver.v[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.dOld[Id(block(x), block(y))] += m_oldFluidSolver.dOld[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.uOld[Id(block(x), block(y))] += m_oldFluidSolver.uOld[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
				m_fluidSolver.vOld[Id(block(x), block(y))] += m_oldFluidSolver.vOld[Id(block(x, oldSize), block(y, oldSize), oldNumBlocksX)];
			}
		}
		for (int x = 0; x < m_numBlocksX; ++x)
		{
			for (int y = 0; y < m_numBlocksY; ++y)
			{
				m_fluidSolver.d[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.u[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.v[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.dOld[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.uOld[Id(x, y)] /= (newSize*newSize);
				m_fluidSolver.vOld[Id(x, y)] /= (newSize*newSize);
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
     * Calculate the mouse input density in the fluid grid. We add density linearly decreasing
     * in circles around the old mouse position (m_MouseXOld, m_MouseYOld).
     *
     * @param pX: density is added ad point with first coordinate pX.
     * @param pY: density is added ad point with first coordinate pY.
     * @param radius: density is added up to radius around
     * old mouse position.
     **/
    private void addDensity(int pX, int pY, int radius)
    {
//    	PsDebug.message("addForce to x: " + String.valueOf(m_mouseX_Old) + ", y: " + String.valueOf(m_mouseY_Old));
//    	PsDebug.message("new x: " + String.valueOf(m_mouseX) + ", new y: " + String.valueOf(m_mouseY));
    	float[] lowerDbound = new float[m_numBlocksX*m_numBlocksY];
    	PdVector mousePos = new PdVector((double)pX, (double)pY);
    	final float densityNormConst = 3/(float)(radius*radius)/(float)Math.PI;
//    	boolean wasZero = true;
    	
    	// Add densities to lower bound, linearly decreasing from center
		for (int x=Math.max(0, pX-radius); x<Math.min(m_imageWidth, pX+radius); x++)
		{
			for (int y=Math.max(0, pY-radius); y<Math.min(m_imageHeight, pY+radius); y++)
			{
//				if (m_fluidSolver.dOld[I(x, y)] != 0)
//					wasZero = false;
//		    	PsDebug.message("dist == " + String.valueOf(mousePos.dist(new PdVector((double)x, (double)y))));
				double lambda = Math.max(0.0, (radius - mousePos.dist(new PdVector((double)x, (double)y))) / (double)radius);
				
//				PsDebug.message("x: " + String.valueOf(x) + ", y: " + String.valueOf(y)
//							+ ", lambda: " + String.valueOf(lambda));

				// m_fluidSolver.dOld[I(x, y)] = Math.max((float)(lambda * m_densityConst.getValue() * densityNormConst), m_fluidSolver.dOld[I(x, y)]);
				lowerDbound[Id(block(x), block(y))] += (float)(lambda * m_densityConst.getValue() * densityNormConst);
			}
		}
		for (int x=Math.max(0, block(pX-radius)); x<Math.min(m_numBlocksX, block(pX+radius)); ++x)
		{
			for (int y=Math.max(0, block(pY-radius)); y<Math.min(m_numBlocksY, block(pY+radius)); ++y)
			{
				// Average
				lowerDbound[Id(x, y)] /= (m_blockSize.getValue() * m_blockSize.getValue());

				// Put maximum of input force that is already there and new forcecone into fluidSolver
				m_fluidSolver.dOld[Id(x, y)] = Math.max(lowerDbound[Id(x, y)], m_fluidSolver.dOld[Id(x, y)]);
			}
		}
		
//		if (!wasZero)
//			PsDebug.warning("dOld war nicht null!");
    }
    
    /**
     * Resets the force- and densityTraces, i.e. deletes every entry except for 
     * the last one, which is then the next "old" mouse position.
     **/
    private void resetMouseTraces()
    {
		// Force traces
//    	int size = m_forceTraceX.getSize();
//    	if (size != m_forceTraceY.getSize())
//    		PsDebug.warning("m_forceTraceX and m_forceTraceX have different sizes!");
//    	else if (size > 0)
//    	{
//			int tempX = m_forceTraceX.getEntry(size-1);
//			int tempY = m_forceTraceY.getEntry(size-1);
//			m_forceTraceX = new PiVector();
//			m_forceTraceX.setSize(1);
//			m_forceTraceX.setEntry(0, tempX);
//			m_forceTraceY = new PiVector();
//			m_forceTraceY.setSize(1);
//			m_forceTraceY.setEntry(0, tempY);
//    	}
		m_forceTraceX = new PiVector();
		m_forceTraceX.setSize(0);
		m_forceTraceY = new PiVector();
		m_forceTraceY.setSize(0);
		
		// Density traces
		m_densityTraceX = new PiVector();
		m_densityTraceX.setSize(0);
		m_densityTraceY = new PiVector();
		m_densityTraceY.setSize(0);
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
	
	// Just used for debugging:
	// DOES NOT YET WORK WITH BLOCKSIZE
	// private void checkWhetherduvOldEqualsZero()
	// {
	// 	// Debugging: Checken, ob uOld, vOld == 0 sind
	// 	for (int x=0; x<m_imageWidth; x+=m_blockSize.getValue())
	// 	{
	// 		boolean breakThis = false;
	// 		for (int y=0; y<m_imageHeight; y+=m_blockSize.getValue())
	// 		{
	// 			if (m_fluidSolver.uOld[I(x,y)] != 0)
	// 			{
	// 				PsDebug.warning("uOld nicht 0 am Beginn eines Frames!");
	// 				breakThis = true;
	// 			}
	// 			else if (m_fluidSolver.vOld[I(x,y)] != 0)
	// 			{
	// 				PsDebug.warning("vOld nicht 0 am Beginn eines Frames!");
	// 				breakThis = true;
	// 			}
	// 			else if (m_fluidSolver.dOld[I(x,y)] != 0)
	// 			{
	// 				PsDebug.warning("dOld nicht 0 am Beginn eines Frames!");
	// 				breakThis = true;
	// 			}
	// 			
	// 			if (breakThis)
	// 				break;
	// 		}
	// 		if (breakThis)
	// 			break;
	// 	}
	// }
	
	public PuDouble getTime()
	{
		return m_time;
	}
}