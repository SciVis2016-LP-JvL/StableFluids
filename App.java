package StableFluids;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.Rectangle;

import jv.object.PsPanel;
import jv.object.PsViewerIf;
import jv.project.PjProject;
import vgp.object.PsApplet;

/**
 * Project visualizing smoke like [Jos Stam 1999].
 * 
 * @author		Lukas Polthier, Johannes von Lindheim, Konrad Polthier (Julia Set project)
 * @version			20.12.15, 0.00 reused JuliaSets project as template for Stable Fluids application
 */
public class App extends PsApplet {
	/** Interface of applet to inform about author, version, and copyright */
	public String getAppletInfo() {
		return "Name: "		+ this.getClass().getName()+ "\r\n" +
				 "Author: "		+ "Konrad Polthier" + "\r\n" +
				 "Version: "	+ "2.00" + "\r\n" +
				 "Applet shows usage of images in display" + "\r\n";
	}
	/** Return a new allocated project instance. */
	public Rectangle getSizeOfFrame() {
//		return new Rectangle(580, 5, 272, 620);
		return new Rectangle(580, 5, 272+256, 720);
	}
	/** Return a new allocated project instance. */
	public PjProject getProject() {
		return new Painter();
	}
	/**
	 * Standalone application support. The main() method acts as the applet's
	 * entry point when it is run as a standalone application. It is ignored
	 * if the applet is run from within an HTML page.
	 */
	public static void main(String args[]) {
		main(new App(), args);
	}
	/**
	 * Configure and initialize the viewer, load system and user projects.
	 * <p>
	 * Overrides the run() method of PsApplet to specify an individual layout.
	 */
	@SuppressWarnings("deprecation")
	public void run() {
		drawMessage("Loading project ...");
		Painter project = (Painter)getProject();
		m_viewer.addProject(project);
		m_viewer.selectProject(project);

		// Must be called after registration of project in m_viewer
		// since m_viewer need in method getDispMandelbrot().
//		project.addDisplay(project.getDispMandelbrot());
		
		// Get 3d display from viewer and add it to applet
		setLayout(new BorderLayout());
		add(project.getInfoPanel(), BorderLayout.CENTER);
		PsPanel pDisplay = new PsPanel(new GridLayout(1, 2));
		{
			pDisplay.setPreferredSize(512, 256);
			pDisplay.add(project.getDisp().getCanvas());
//			pDisplay.add(project.getDispMandelbrot().getCanvas());
		}
		add(pDisplay, BorderLayout.NORTH);
		validate();

		// Choose initial panel in control window (must press F1 inside the applet)
		m_viewer.showPanel(PsViewerIf.MATERIAL);

		// Explicitly start the applet
		startFromThread();
	}
}
