package StableFluids;

//import java.awt.BorderLayout;
import java.awt.Button;
//import java.awt.Canvas;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Panel;
import java.awt.Label;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

//import jv.number.PdVector_IP;
//import jv.object.PsConfig;
//import jv.object.PsMainFrame;
import jv.object.PsPanel;
import jv.object.PsUpdateIf;
import jv.objectGui.PsMultiLineLabel;
import jv.project.PjProject_IP;
//import jv.project.PvDisplayIf;
//import jv.vecmath.PdVector;

/**
 * Project visualizing smoke like [Jos Stam 1999].
 * 
 * @author		Lukas Polthier, Johannes von Lindheim, [Konrad Polthier (Julia Set project)]
 * @version			20.12.15, 0.00 reused JuliaSets project as template for Stable Fluids application
 */
public class Painter_IP extends PjProject_IP implements ActionListener { 
	protected	Painter				m_Project;
	protected	PsPanel				m_pSlider;
	protected	Panel				m_pBottomButtons;
	protected	Button				m_bReset;
	protected	Button				m_bClear;
	protected 	Button				m_bFlipColor;
	protected 	Button				m_bFreeze;
	protected 	Button				m_bImportImage;
	protected	PsPanel				m_lTime;
	
	protected final String notice	= "Drag mouse:\n"
			+ "Left - insert smoke\n"
			+ "Right - add velocity";

	public Painter_IP() {
		super();

		if (getClass() == Painter_IP.class) {
			init();
		}
	}
	public void init() {
		super.init();
		addTitle("");

		// Description
		PsPanel pNotice = new PsPanel();
		pNotice.setInsetSize(5);
		pNotice.setBorderType(PsPanel.BORDER_GROOVE);
		pNotice.add(new PsMultiLineLabel(notice));
		add(pNotice);
		
		// Time
// 		m_Project.getTime().addInspector("Time Inspector", m_lTime);
// 		Panel pTime = new Panel();
// 		pTime.setLayout(new GridLayout(1, 2));
// 		pTime.add(new Label("Time"));
// 		m_lTime = new PsPanel();
// 		pTime.add(m_lTime);
		
		// Slider Panel
		m_pSlider = new PsPanel();
		add(m_pSlider);
		
		// buttons at bottom
		m_pBottomButtons = new Panel();
		m_pBottomButtons.setLayout(new FlowLayout(FlowLayout.CENTER));
		add(m_pBottomButtons);
		m_bReset = new Button("Reset App");
		m_bReset.addActionListener(this);
		m_pBottomButtons.add(m_bReset);
		
		m_bClear = new Button("Clear Canvas");
		m_bClear.addActionListener(this);
		m_pBottomButtons.add(m_bClear);
		
		m_bFlipColor = new Button("Flip Color");
		m_bFlipColor.addActionListener(this);
		m_pBottomButtons.add(m_bFlipColor);
		
		m_bImportImage = new Button("Import Image");
		m_bImportImage.addActionListener(this);
		m_pBottomButtons.add(m_bImportImage);

		m_bFreeze = new Button("Freeze");
		m_bFreeze.addActionListener(this);
		m_pBottomButtons.add(m_bFreeze);
	}
	
	/**
	 * Set parent of panel which supplies the data inspected by the panel.
	 */
	public void setParent(PsUpdateIf parent) {
		super.setParent(parent);
		m_Project = (Painter)parent;
		setTitle("Stable Fluids");
		m_pSlider.add(m_Project.m_blockSize.getInfoPanel());
		m_pSlider.add(m_Project.m_densityRadius.getInfoPanel());
		m_pSlider.add(m_Project.m_densityConst.getInfoPanel());
		m_pSlider.add(m_Project.m_forceRadius.getInfoPanel());
		m_pSlider.add(m_Project.m_forceConst.getInfoPanel());
		m_pSlider.add(m_Project.m_buoyancy.getInfoPanel());
		m_pSlider.add(m_Project.m_diffusion.getInfoPanel());
		m_pSlider.add(m_Project.m_viscosity.getInfoPanel());
		m_pSlider.add(m_Project.m_vorticity.getInfoPanel());
	}
	
	/**
	 * Update the panel whenever the parent has changed somewhere else.
	 * Method is invoked from the parent or its superclasses.
	 */
	public boolean update(Object event) {
		if (m_Project == event) {
			return true;
		}
		return super.update(event);
	}
	
	/**
	 * Handle action events invoked from buttons, menu items, text fields.
	 */
	public void actionPerformed(ActionEvent event) {
		if (m_Project==null)
			return;
		Object source = event.getSource();
		if (source == m_bReset) {
//			m_Project = new Painter();
			// m_Project.reset();
			// m_Project.init();
			// m_Project.start();
			m_Project.buttonReset();
			// m_pSlider.update(this);
			m_Project.update(m_Project);
		}
		if (source == m_bClear) {
			m_Project.buttonClear();
			m_Project.update(m_Project);
		}
		if (source == m_bFlipColor) {
			m_Project.buttonFlipColor();
			m_Project.update(m_Project);
		}
		if (source == m_bImportImage) {
			m_Project.buttonImportImage();
			m_Project.update(m_Project);
		}
		if (source == m_bFreeze) {
			m_Project.buttonFreeze();
			m_Project.update(m_Project);
		}
	}
}

