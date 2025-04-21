import tkinter as tk
from tkinter import filedialog, Frame, Label, Button, Scale, DoubleVar, HORIZONTAL, Checkbutton
import pandas as pd
import numpy as np
from colorsys import rgb_to_hls
from scipy.interpolate import griddata
from scipy.optimize import minimize, curve_fit
import matplotlib.colors as mcolors
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import messagebox
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import griddata
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
from sklearn.metrics import r2_score  

mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Times New Roman'

# Global variables 
root = tk.Tk()
# root.title("Frequency-Temperature Superposition Tool")
root.title("Temporal-Modulus Analyzer v1.1.4")
# root.geometry("1400x950")
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()
init_width = int(screen_width * 0.82)
init_height = int(screen_height * 0.82)
root.geometry(f"{init_width}x{init_height}")
root.minsize(1400, 720)

data = None
analysis_shift_factors = {}   # Shift factor calculated on page 3
fitted_params = {}            # Page 4.Fitted parameters {'a':... , 'b':... , 'c':... , 'd':...  }
abcd_parameters_global = None

# Global canvas for each page (if you need to destroy old canvases)
canvas_raw = None       # Page 1: Raw data plot
canvas_3d = None        # Page 2:3D Surface plot
canvas_analysis = None  # Page 3: Analysis diagram
canvas_curve_fit = None # Page 4: Fit plot
canvas_eprime = None    # Page 5: E'(w) graph
canvas_et = None        # Page 6: E(t) plot
canvas_em = None        # Page 8: Elastic modulus vs strain rate plot
canvas_mod_temp = None  # Page 11: Modulus vs Temperature
current_fig = None

# Page 3
scale_var = None
analysis_frame = None
shift_factors_label = None

# Page 4
var_a_upper = None
var_d_upper = None
fit_params_label = None

# Global strain parameter
strain_min = 1e-25
strain_max = 0.0025
num_steps = 500
strain_rates_to_plot = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]


###############################################
# Help function
###############################################
def is_dark_color(hex_color):
    rgb = mcolors.hex2color(hex_color)
    h, l, s = rgb_to_hls(*rgb)
    return l < 0.5

def add_watermark(ax, text="NYU-ViscoMOD", on_axes=False):
    if on_axes:        
        ax.figure.text(
            0.99, 0.01,  # x=0.99, y=0.01 
            text,
            fontsize=14,
            fontname="Times New Roman",
            color="purple",
            ha="right", va="bottom", alpha=0.5,
            transform=ax.figure.transFigure  
        )
    else:        
        ax.text(
            0.99, 0.01,  
            text,
            fontsize=14,
            fontname="Times New Roman",
            color="purple",
            ha="right", va="bottom", alpha=0.5,
            transform=ax.transAxes  
        )


###############################################
# Page 0: Load CSV 
###############################################
tension_data = None
torsion_data = None

def show_about():
    message = (
        "Copyright to NYU\n"
        "Patent number US Patent #10,345,210\n"
        "The software is freely available for use and publication, "
        "provided that appropriate acknowledgement is given.\n"
        "*Disclaimer: The developer takes no responsibility of the calculation."
    )
    messagebox.showinfo("About", message)

def load_csv_page(content_frame):
    global tension_data, torsion_data 
    tension_data = None
    torsion_data = None

    for widget in content_frame.winfo_children():
        widget.destroy()
    
    page_frame = Frame(content_frame)
    page_frame.pack(fill="both", expand=True)

    toolbar = Frame(page_frame)
    toolbar.pack(side="top", anchor="nw", fill="x")
    help_mb = tk.Menubutton(toolbar, text="Help", font=('Times New Roman', 14))
    help_menu = tk.Menu(help_mb, tearoff=0)
    help_menu.add_command(label="About", command=show_about, font=('Times New Roman', 12))
    help_mb.config(menu=help_menu)
    help_mb.pack(side="left", padx=5, pady=5)

    left_frame = Frame(page_frame)
    left_frame.pack(padx=5, pady=5, fill='both', expand=True, anchor='center')

    left_frame.grid_columnconfigure(0, weight=1)
    left_frame.grid_columnconfigure(1, weight=1)

    lbl = Label(left_frame, text="Select CSV files to load data:", font=('Times New Roman', 22))
    lbl.grid(row=0, column=0, columnspan=2, pady=5)

    # Tension0
    def load_tension0_csv():
        global tension_data
        fp = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if fp:
            try:
                tension_data = pd.read_csv(fp)
                tension_data['Frequency'] = tension_data['Frequency'].astype(float)
                tension_data['Storage Modulus'] = tension_data['Storage Modulus'].astype(float)
                messagebox.showinfo("Info", "Tension data loaded successfully.")
                set_global_data()
            except Exception as e:
                messagebox.showerror("Error", "Load failure: " + str(e)) 

    # Tension
    def load_tension_csv():
        global tension_data
        fp = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if fp:
            try:
                tension_data = pd.read_csv(fp)
                tension_data['Frequency'] = tension_data['Frequency'].astype(float)
                tension_data['Storage Modulus'] = tension_data['Storage Modulus'].astype(float)
                messagebox.showinfo("Info", "Tension data loaded successfully.")
                set_global_data()
            except Exception as e:
                messagebox.showerror("Error", "Load failure: " + str(e))
    
    # Torsion
    def load_torsion_csv():
        global torsion_data
        fp = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if fp:
            try:
                torsion_data = pd.read_csv(fp)
                torsion_data['Frequency'] = torsion_data['Frequency'].astype(float)
                torsion_data['Storage Modulus'] = torsion_data['Storage Modulus'].astype(float)
                messagebox.showinfo("Info", "Torsion data loaded successfully.")
                set_global_data()
            except Exception as e:
                messagebox.showerror("Error", "Load failure: " + str(e))

    # Tension1
    def load_tension1_csv():
        global tension_data
        fp = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if fp:
            try:
                tension_data = pd.read_csv(fp)
                tension_data['Frequency'] = tension_data['Frequency'].astype(float)
                tension_data['Storage Modulus'] = tension_data['Storage Modulus'].astype(float)
                messagebox.showinfo("Info", "Tension data loaded successfully.")
                set_global_data()
            except Exception as e:
                messagebox.showerror("Error", "Load failure: " + str(e)) 
        
    btn_load_tension = Button(left_frame, text="Load Tension CSV", font=('Times New Roman', 22), command=load_tension_csv, width=30, height=5)
    btn_load_tension.grid(row=1, column=1, padx=(50,20), pady=20, sticky='nsew')
    btn_load_tension0 = Button(left_frame, text="Load Cantilever CSV", font=('Times New Roman', 22), command=load_tension_csv, width=30, height=5)
    btn_load_tension0.grid(row=1, column=0, padx=(50,20), pady=20, sticky='ew')

    btn_load_tension2 = Button(left_frame, text="Load Compression CSV", font=('Times New Roman', 22), command=load_tension_csv, width=30, height=5)
    btn_load_tension2.grid(row=2, column=0, padx=(50,20), pady=20, sticky='ew')
    btn_load_torsion = Button(left_frame, text="Load Torsion CSV", font=('Times New Roman', 22), command=load_torsion_csv, width=30, height=5)
    btn_load_torsion.grid(row=2, column=1, padx=(50,20), pady=20, sticky='ew')

    btn_next = Button(left_frame, text="Next", font=('Times New Roman', 22), command=lambda: show_page(raw_data_page))
    btn_next.grid(row=3, column=0, columnspan=2, pady=20)

    # Added label at the bottom for user manual note
    lbl_manual = Label(left_frame, text="*Please read the user manual before using the software", font=('Times New Roman', 20))
    lbl_manual.grid(row=4, column=0, columnspan=2, pady=20)

def set_global_data():
    global tension_data, torsion_data, data
    if tension_data is not None:
        data = tension_data
    elif torsion_data is not None:
        data = torsion_data
    else:
        data = None


###############################################
# Page 1: Raw Data
###############################################
def raw_data_page(content_frame):
 
    for widget in content_frame.winfo_children():
        widget.destroy()

    global data, canvas_raw
    if data is None:
        Label(content_frame, text="No data loaded!", font=('Times New Roman', 12)).pack()
        return

    top_frame = Frame(content_frame)
    top_frame.pack(side="top", fill="x", padx=10, pady=10)

    btn_prev = Button(top_frame, text="Previous", font=('Times New Roman', 12),
                      command=lambda: show_page(load_csv_page))
    btn_prev.pack(side="left", padx=5)


    btn_next = Button(top_frame, text="Next", font=('Times New Roman', 12),
                      command=lambda: show_page(show_surface_plot))
    btn_next.pack(side="right", padx=5)

    btn_save_graph = Button(top_frame, text="Save Graph", font=('Times New Roman', 12),
                            command=lambda: save_graph(canvas_raw))
    btn_save_graph.pack(side="right", padx=5)

    main_frame = Frame(content_frame)
    main_frame.pack(side="top", fill="both", expand=True)

    fig = Figure(figsize=(8, 5))
    ax = fig.add_subplot(111)

    temps = data['Temperature'].unique()
    all_colors = list(mcolors.CSS4_COLORS.values())
    darker = [c for c in all_colors if is_dark_color(c)]

    for i, temp in enumerate(temps):
        subdata = data[data['Temperature'] == temp]
        ax.semilogx(subdata['Frequency'], subdata['Storage Modulus'],
                    color=darker[i % len(darker)], label=f"{temp} °C")

    ax.set_xscale('log')
    xmin_val = data['Frequency'].min()
    xmax_val = data['Frequency'].max()
    ax.set_xlim(xmin_val * 0.8, xmax_val * 1.2)

    ymin_val = data['Storage Modulus'].min()
    ymax_val = max(data['Storage Modulus'].max(), 10)
    lower_bound = 0 if ymin_val < 0 else ymin_val * 0
    ax.set_ylim(lower_bound, ymax_val * 1.1)

    ax.set_xlabel('Frequency (Hz)', fontname="Times New Roman", fontsize=20)
    ax.set_ylabel("Storage Modulus (MPa)", fontname="Times New Roman", fontsize=20)
    ax.set_title('Raw Data Plot', fontname="Times New Roman", fontsize=20)
    fig.tight_layout()
    fig.subplots_adjust(right=0.8)

    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontname("Times New Roman")
        label.set_fontsize(16)

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),
              prop={"family": "Times New Roman", "size": 12})

    add_watermark(ax, "NYU-ViscoMOD", on_axes=True)

    if canvas_raw:
        canvas_raw.get_tk_widget().destroy()

    canvas_raw = FigureCanvasTkAgg(fig, master=main_frame)
    canvas_raw.draw()
    canvas_raw.get_tk_widget().pack(fill="both", expand=True)

# Save graph function
def save_graph(canvas):
    fig = canvas.figure 
    file_path = filedialog.asksaveasfilename(defaultextension=".jpeg", filetypes=[("JPEG files", "*.jpeg")])
    if file_path:

        fig.savefig(file_path, dpi=300, format='jpeg', bbox_inches='tight') 
        messagebox.showinfo("Info", f"Graph saved to {file_path}")


###############################################
# Page 2: 3D Surface
###############################################
def show_surface_plot(content_frame):
    if data is None or data.empty:
        messagebox.showerror("Error","No data loaded!")
        return
    for widget in content_frame.winfo_children():
        widget.destroy()
    page_frame = Frame(content_frame)
    page_frame.pack(fill="both", expand=True)
    top_frame = Frame(page_frame)
    top_frame.pack(side="top", fill="x", padx=5, pady=5)
    btn_back = Button(top_frame, text="Previous", font=('Times New Roman', 12), command=lambda: show_page(raw_data_page))
    btn_back.pack(side="left", padx=5)

    right_frame = Frame(top_frame)
    right_frame.pack(side="right")
    
    # Save Graph button
    btn_save_graph = Button(right_frame, text="Save Graph", font=('Times New Roman', 12), command=lambda: save_graph(canvas_3d))
    btn_save_graph.pack(side="left", padx=5)
    
    # Next button
    btn_next = Button(right_frame, text="Next", font=('Times New Roman', 12), command=lambda: show_page(analysis_page))
    btn_next.pack(side="left", padx=5)
    canvas_frame = Frame(page_frame)
    canvas_frame.pack(fill="both", expand=True)
    X = data['Temperature']
    Y = data['Frequency']
    Z = data['Storage Modulus']
    x_grid, y_grid = np.meshgrid(np.linspace(X.max(), X.min(), 100),
                                 np.linspace(Y.min(), Y.max(), 100))
    Z_grid = griddata((X, Y), Z, (x_grid, y_grid), method='cubic')
    fig = Figure(figsize=(8,5))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(x_grid, y_grid, Z_grid, cmap='rainbow', edgecolor='none')
    ax.set_xlabel('Temperature (°C)', fontname="Times New Roman", fontsize=18)
    ax.set_ylabel('Frequency (Hz)', fontname="Times New Roman", fontsize=18)
    ax.set_zlabel('Storage Modulus (MPa)', fontname="Times New Roman", fontsize=18)

    for label in ax.get_xticklabels() + ax.get_yticklabels() + ax.zaxis.get_ticklabels():
        label.set_fontname("Times New Roman")
        label.set_fontsize(14)
    # cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.08)
    # add colorbar
    cax = fig.add_axes([0.86, 0.15, 0.03, 0.7])  # (left, bottom, width, height)
    color_bar = fig.colorbar(surf, cax=cax, aspect=20)
    color_bar.set_label('Storage Modulus (MPa)', fontname="Times New Roman", fontsize=20)
    for t in color_bar.ax.get_yticklabels():
        t.set_fontname('Times New Roman')
        t.set_fontsize(14)

    fig.subplots_adjust(left=0.08, right=0.9, bottom=0.1, top=0.9)
    add_watermark(ax, "NYU-ViscoMOD", on_axes=True)
    canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill="both", expand=True)
    global canvas_3d
    canvas_3d = canvas
    # add_watermark(content_frame, "NYU-ViscoMOD")


# ###############################################
# # Page 3: Analysis (live lower_Hz slider)
# ###############################################
import tkinter as tk
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import pandas as pd
import numpy as np
from scipy.optimize import minimize

analysis_shift_factors = {}

def analysis_page(content_frame):
    if data is None or data.empty:
        messagebox.showerror("Error", "No data loaded!")
        return

    for widget in content_frame.winfo_children():
        widget.destroy()

    page_frame = tk.Frame(content_frame)
    page_frame.pack(fill="both", expand=True)

    top_frame = tk.Frame(page_frame)
    top_frame.pack(side="top", fill="x", padx=5, pady=5)

    btn_back = tk.Button(top_frame, text="Previous", font=('Times New Roman', 12), command=lambda: show_page(show_surface_plot))
    btn_back.pack(side="left", padx=5)

    right_frame = tk.Frame(top_frame)
    right_frame.pack(side="right")
    
    btn_save_graph = tk.Button(right_frame, text="Save Graph", font=('Times New Roman', 12), command=lambda: save_graph(canvas_analysis))
    btn_save_graph.pack(side="left", padx=5)
    
    # btn_next = tk.Button(right_frame, text="Next", font=('Times New Roman', 12), command=lambda: show_page(curve_fitting_page))
    btn_next = tk.Button(right_frame, text="Next", font=('Times New Roman', 12), command=go_to_curve_fitting_page)
    btn_next.pack(side="left", padx=5)

    analysis_frame_local = tk.Frame(page_frame)
    analysis_frame_local.pack(fill="both", expand=True)
    lbl_sf = tk.Label(analysis_frame_local, text="", justify="left", font=('Times New Roman', 10))
    lbl_sf.pack(side="right", padx=5)
    
    global analysis_frame, shift_factors_label
    analysis_frame = analysis_frame_local
    shift_factors_label = lbl_sf

    on_slider_change(content_frame)

def on_slider_change(content_frame):
    global analysis_shift_factors

    if data is None or data.empty:
        return

    reference_temp = min(data["Temperature"].unique())
    df_ref = data[data["Temperature"] == reference_temp]

    extended_data = df_ref.copy()

    temperatures = sorted(data["Temperature"].unique())
    temperatures = [temp for temp in temperatures if temp > reference_temp]

    shift_factors = {reference_temp: 1.0}

    for temp in temperatures:
        df_temp = data[data["Temperature"] == temp]

        max_freq = df_temp["Frequency"].max()
        modulus_max_freq = df_temp.loc[df_temp["Frequency"] == max_freq, "Storage Modulus"].values[0]

        closest_match = extended_data.iloc[(extended_data["Storage Modulus"] - modulus_max_freq).abs().argmin()]
        shift_factor = closest_match["Frequency"] / max_freq  

        shift_factors[temp] = shift_factor

        shifted_freq = df_temp["Frequency"] * shift_factor
        shifted_modulus = df_temp["Storage Modulus"]

        shifted_df = pd.DataFrame({"Frequency": shifted_freq, "Storage Modulus": shifted_modulus, "Temperature": temp})
        extended_data = pd.concat([extended_data, shifted_df], ignore_index=True)

    fig = Figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    fig.subplots_adjust(right=0.8)

    ax.plot(df_ref["Frequency"], df_ref["Storage Modulus"], '--', label=f"{reference_temp}°C Original", linewidth=2)

    for temp in temperatures:
        shifted_df = extended_data[extended_data["Temperature"] == temp]
        ax.plot(shifted_df["Frequency"], shifted_df["Storage Modulus"], '--', label=f"{temp}°C Shifted", linewidth=2)

    ax.set_xscale("log")
    ax.set_yscale("linear")
    ax.set_xlabel("Frequency (Hz)", fontname="Times New Roman", fontsize=20)
    ax.set_ylabel("Storage Modulus (MPa)", fontname="Times New Roman", fontsize=20)
    # ax.set_title(f"Extended {reference_temp}°C Data Using Horizontally Shifted Higher Temperature Data", fontname="Times New Roman", fontsize=20)
    ax.set_title(f"Master Curve at {reference_temp}°C", fontname="Times New Roman", fontsize=20)
    ax.grid(True, which="both", linestyle="--", linewidth=0.5)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontname("Times New Roman")
        label.set_fontsize(16)
    
    ymin_val = data['Storage Modulus'].min()
    ymax_val = max(data['Storage Modulus'].max(), 10)
    if ymin_val < 0:
        lower_bound = 0
    else:
        lower_bound = ymin_val * 0
    ax.set_ylim(lower_bound, ymax_val * 1.1)

    # fig.tight_layout()
    add_watermark(ax, "NYU-ViscoMOD", on_axes=True)
    
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), prop={"family": "Times New Roman", "size": 12})

    global canvas_analysis
    if canvas_analysis:
        canvas_analysis.get_tk_widget().destroy()
    canvas_analysis = FigureCanvasTkAgg(fig, master=analysis_frame)
    canvas_analysis.draw()
    canvas_analysis.get_tk_widget().pack(side="left", fill="both", expand=True)

    shift_lines = [f"Shift Factors:"]
    for t, s in shift_factors.items():
        shift_lines.append(f"  {t}°C : {s:.6e}")
    shift_factors_label.config(text="\n".join(shift_lines), font=("Times New Roman", 14))
    analysis_shift_factors.clear()
    analysis_shift_factors.update(shift_factors)

def show_page(page_func):
    content_frame = tk.Frame(root)
    page_func(content_frame)


###################################################
# Page 4: Curve Fitting (Upper Bound for 'a' & 'd')
###################################################
def update_slider_page4(var, delta, callback, var_a, var_d):

    val = var.get()
    new_val = min(max(val + delta, 0), 2000)  
    var.set(new_val)

    callback(var_a, var_d)

def go_to_analysis_page():
    global canvas_curve_fit
    canvas_curve_fit = None
    show_page(analysis_page)

def go_to_master_curve_each_temp_page():
    global canvas_curve_fit
    canvas_curve_fit = None
    show_page(master_curve_each_temp_page)

def go_to_curve_fitting_page():
    global canvas_curve_fit
    canvas_curve_fit = None
    show_page(curve_fitting_page)

def curve_fitting_page(content_frame, shift_factors=None):
    global analysis_shift_factors, var_a_upper_local, var_d_upper_local, df_results_global, df_results_global_page12, abcd_parameters_global
    df_results_global = None
    df_results_global_page12 = None
    abcd_parameters_global = []
    if data is None or data.empty:
        messagebox.showerror("Error","No data loaded!")
        return
    if not analysis_shift_factors:
        messagebox.showerror("Error","No shift factors found. Please complete Page 3 first.")
        return
    for widget in content_frame.winfo_children():
        widget.destroy()
    page_frame = Frame(content_frame)
    page_frame.pack(fill="both", expand=True)
    top_frame = Frame(page_frame)
    top_frame.pack(side="top", fill="x", padx=5, pady=5)
    # btn_back = Button(top_frame, text="Previous", font=('Times New Roman', 12), command=lambda: show_page(analysis_page))
    btn_back = Button(top_frame, text="Previous", font=('Times New Roman', 12), command=go_to_analysis_page)
    btn_back.pack(side="left", padx=5)
    lbl_a = Label(top_frame, text="Upper Bound for 'a':", font=('Times New Roman', 12))
    lbl_a.pack(side="left", padx=5)
    var_a_upper_local = DoubleVar(value=500.0)
    s_a = Scale(top_frame, from_=0.0, to=2000.0, resolution=1.0,
                orient=HORIZONTAL, length=200,
                variable=var_a_upper_local,
                command=lambda v: on_curve_fit_bounds_changed(var_a_upper_local, var_d_upper_local))
    s_a.pack(side="left", padx=2)
    btn_a_minus = Button(top_frame, text="-", font=('Times New Roman', 12), command=lambda: update_slider_page4(var_a_upper_local, -1.0, on_curve_fit_bounds_changed, var_a_upper_local, var_d_upper_local))
    btn_a_minus.pack(side="left", padx=2)
    btn_a_plus = Button(top_frame, text="+", font=('Times New Roman', 12), command=lambda: update_slider_page4(var_a_upper_local, 1.0, on_curve_fit_bounds_changed, var_a_upper_local, var_d_upper_local))
    btn_a_plus.pack(side="left", padx=2)
    lbl_d = Label(top_frame, text="  Upper Bound for 'd':", font=('Times New Roman', 12))
    lbl_d.pack(side="left", padx=5)
    var_d_upper_local = DoubleVar(value=500.0)
    s_d = Scale(top_frame, from_=0.0, to=2000.0, resolution=1.0,
                orient=HORIZONTAL, length=200,
                variable=var_d_upper_local,
                command=lambda v: on_curve_fit_bounds_changed(var_a_upper_local, var_d_upper_local))
    s_d.pack(side="left", padx=2)
    btn_d_minus = Button(top_frame, text="-", font=('Times New Roman', 12), command=lambda: update_slider_page4(var_d_upper_local, -1.0, on_curve_fit_bounds_changed, var_a_upper_local, var_d_upper_local))
    btn_d_minus.pack(side="left", padx=2)
    btn_d_plus = Button(top_frame, text="+", font=('Times New Roman', 12), command=lambda: update_slider_page4(var_d_upper_local, 1.0, on_curve_fit_bounds_changed, var_a_upper_local, var_d_upper_local))
    btn_d_plus.pack(side="left", padx=2)
    right_frame = Frame(top_frame)
    right_frame.pack(side="right")
    btn_save_graph = Button(right_frame, text="Save Graph", font=('Times New Roman', 12), command=lambda: save_graph(canvas_curve_fit))
    btn_save_graph.pack(side="left", padx=5)
    # btn_next = Button(right_frame, text="Next", font=('Times New Roman', 12), command=lambda: show_page(master_curve_each_temp_page))
    btn_next = Button(right_frame, text="Next", font=('Times New Roman', 12), command=go_to_master_curve_each_temp_page)
    btn_next.pack(side="left", padx=5)
    # btn_next = Button(top_frame, text="Next", font=('Times New Roman', 12), command=lambda: show_page(eprime_plot_page))
    # btn_next.pack(side="right", padx=5)
    curve_frame = Frame(page_frame)
    curve_frame.pack(fill="both", expand=True)
    lbl_fit = Label(curve_frame, text="E'(w)= a tanh(b*(w+c)) + d", justify="left", font=('Times New Roman', 14))
    lbl_fit.pack(side="right", padx=10)
    global fit_params_label
    fit_params_label = lbl_fit
    on_curve_fit_bounds_changed(var_a_upper_local, var_d_upper_local)

def on_curve_fit_bounds_changed(var_a_upper_local, var_d_upper_local, extra=None):
    do_curve_fitting(var_a_upper_local.get(), var_d_upper_local.get())

def do_curve_fitting(a_upper_val, d_upper_val):
    global data, analysis_shift_factors, canvas_curve_fit, fit_params_label, fitted_params
    if data is None or data.empty:
        return
    if not analysis_shift_factors:
        return
    combined_log_freq = []
    combined_modulus = []
    temps_sorted = sorted(data['Temperature'].unique())
    for t in temps_sorted:
        subdata = data[data['Temperature'] == t]
        sf = analysis_shift_factors.get(t, 1.0)
        shifted_log_f = np.log10(subdata['Frequency'] * sf)
        combined_log_freq.extend(shifted_log_f)
        combined_modulus.extend(subdata['Storage Modulus'])
    combined_log_freq = np.array(combined_log_freq)
    combined_modulus = np.array(combined_modulus)

    def storage_modulus_model(log_omega, a, b, c, d):
        return a * np.tanh(b*(log_omega + c)) + d
    lower_bounds = [1e-6, -100.0, -100.0, 1e-6]
    # upper_bounds = [a_upper_val, 100.0, 100.0, d_upper_val]
    upper_bounds = [a_upper_val, 100.0, 100.0, d_upper_val]
    try:
        params, _ = curve_fit(storage_modulus_model,
                              combined_log_freq, combined_modulus,
                              bounds=(lower_bounds, upper_bounds),
                              maxfev=1000000)
    except Exception as e:
        print("Curve fitting error:", e)
        return
    a, b, c, d = params
    fitted_params['a'] = a
    fitted_params['b'] = b
    fitted_params['c'] = c
    fitted_params['d'] = d
    # fig = Figure(figsize=(8,5))

    if canvas_curve_fit is None:
        fig = Figure(figsize=(8, 5))
        canvas_curve_fit = FigureCanvasTkAgg(fig, master=fit_params_label.master)
        canvas_curve_fit.get_tk_widget().pack(side="left", fill="both", expand=True)
    else:
        fig = canvas_curve_fit.figure
        fig.clf() 
    ax = fig.add_subplot(111)
    all_colors = list(mcolors.CSS4_COLORS.values())
    dark_colors = [clr for clr in all_colors if is_dark_color(clr)]
    for i, t in enumerate(temps_sorted):
        sf = analysis_shift_factors.get(t, 1.0)
        subdata = data[data['Temperature'] == t]
        logx = np.log10(subdata['Frequency'] * sf)
        color = dark_colors[i % len(dark_colors)]
        ax.scatter(logx, subdata['Storage Modulus'], color=color, s=10,
                   label="Master Curve" if i==0 else None)
    logmin, logmax = np.min(combined_log_freq), np.max(combined_log_freq)
    log_axis = np.linspace(logmin, logmax, 1000)
    yfit = storage_modulus_model(log_axis, a, b, c, d)
    ax.plot(log_axis, yfit, 'r-', linewidth=2, label="Fitted Master Curve")
    ax.set_xlabel("Reduced Frequency Log(Hz)", fontname="Times New Roman", fontsize=20)
    ax.set_ylabel("Storage Modulus (MPa)", fontname="Times New Roman", fontsize=20)
    ax.set_title("Master Curve with Fitted Equation", fontname="Times New Roman", fontsize=20)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontname("Times New Roman")
        label.set_fontsize(16)

    ymin_val = np.min(combined_modulus)
    ymax_val = np.max(combined_modulus)
    if ymin_val < 0:
        lower_bound = 0
    else:
        lower_bound = ymin_val * 0
    ax.set_ylim(lower_bound, ymax_val * 1.1)

    ax.legend(fontsize=12, prop={"family": "Times New Roman"})
    ax.grid(True)
    add_watermark(ax, "NYU-ViscoMOD", on_axes=True)
    fig.tight_layout()

    yfit = storage_modulus_model(combined_log_freq, a, b, c, d)
    r_squared = r2_score(combined_modulus, yfit)

    canvas_curve_fit.draw()
    txt = (f"E'(w)= a tanh(b*(w+c)) + d\n"
           f"a = {a:.6f}\n"
           f"b = {b:.6f}\n"
           f"c = {c:.6f}\n"
           f"d = {d:.6f}\n"
           f"(Bounds: a ≤ {a_upper_val:.4g}, d ≤ {d_upper_val:.4g})\n"
           f"\nR² = {r_squared:.6f}\n")
    if fit_params_label:
        fit_params_label.config(text=txt)


###############################################
# Define Etime_time_cycle (corrected version)
###############################################
def Etime_time_cycle(time, cycle, a, b, c, d):
    Etime = np.zeros_like(time)
    N1, N2, N3 = 240, 74, 24
    def E_prime(w, a, b, c, d):
        return a * np.tanh(b * (np.log(w) + c)) + d   
    def integrand(t_val, E_prime_w, w_val):
        return (2/np.pi) * (E_prime_w / w_val) * np.sin(w_val * t_val)
    for i, t_val in enumerate(time):
        w1 = np.linspace(1e-6 / t_val, cycle * 0.1 * 2 * np.pi / t_val, int(cycle * 0.1 * N1) + 1)
        w2 = np.linspace(cycle * 0.1 * 2 * np.pi / t_val, cycle * 0.4 * 2 * np.pi / t_val, int(cycle * 0.3 * N2) + 1)
        w3 = np.linspace(cycle * 0.4 * 2 * np.pi / t_val, cycle * 2 * np.pi / t_val, int(cycle * 0.6 * N3) + 1)
        all_w = np.concatenate([w1, w2[1:], w3[1:]])
        y = integrand(t_val, E_prime(all_w, a, b, c, d), all_w)
        Etime[i] = np.trapz(y, all_w)
    return Etime


####################################################################
# Page 5: Master Curve for each Reference Temp (Fitting Parameters)
####################################################################
import numpy as np
import matplotlib.cm as cm
from tkinter import Scale, Button, Frame, Tk, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from matplotlib.figure import Figure

def storage_modulus_model(log_omega, a, b, c, d):
    return a * np.tanh(b * (log_omega + c)) + d

def update_slider(var, delta, callback):
    val = var.get()
    new_val = min(max(val + delta, 0), 2000)  
    var.set(new_val)  
    callback() 

def master_curve_each_temp_page(content_frame):
    global abcd_parameters_global, df_results_global_page12, df_results_global
    df_results_global = None
    df_results_global_page12 = None
    if data is None or data.empty:
        messagebox.showerror("Error", "No data loaded!")
        return
    if not analysis_shift_factors:
        messagebox.showerror("Error", "No shift factors available. Please complete Page 3 first.")
        return

    original_shift_factors = analysis_shift_factors
    estimated_shift_factors = {}
    for ref_temp in original_shift_factors.keys():
        scale_factor = 1 / original_shift_factors[ref_temp]
        estimated_shift_factors[ref_temp] = {temp: factor * scale_factor for temp, factor in original_shift_factors.items()}
    temperatures = np.sort(data['Temperature'].unique())

    abcd_parameters_global = [] 

    for widget in content_frame.winfo_children():
        widget.destroy()
    page_frame = Frame(content_frame)
    page_frame.pack(fill="both", expand=True)
    top_frame = Frame(page_frame)
    top_frame.pack(side="top", fill="x", padx=5, pady=5)

    # btn_prev = Button(top_frame, text="Previous", font=('Times New Roman', 12), command=lambda: show_page(curve_fitting_page))
    btn_prev = Button(top_frame, text="Previous", font=('Times New Roman', 12), command=go_to_curve_fitting_page)
    btn_prev.pack(side="left", padx=5)
    right_frame = Frame(top_frame)
    right_frame.pack(side="right")
    btn_save_graph = Button(right_frame, text="Save Graph", font=('Times New Roman', 12), command=lambda: save_graph(canvas_10))
    btn_save_graph.pack(side="left", padx=5)
    btn_next = Button(right_frame, text="Next", font=('Times New Roman', 12), command=lambda: show_page(elastic_modulus_vs_strain_rate_page))
    btn_next.pack(side="left", padx=5)

    fig = Figure(figsize=(12, 8))
    ax = fig.add_subplot(111)

    n_ref_temps = len(estimated_shift_factors.keys())
    colors = cm.rainbow(np.linspace(0, 1, n_ref_temps))

    a_fit_1 = var_a_upper_local.get()
    d_fit_1 = var_d_upper_local.get()

    b_fit, c_fit = 1, 0

    def update_graph(a, d):
        global abcd_parameters_global
        abcd_parameters_global = []  
        abcd_parameters = []  
        a_fit = a
        d_fit = d
        ax.clear()
        for t in ax.texts[:]:
            t.remove()
        for t in fig.texts[:]:
            t.remove()

        fit_line_added = False

        for ref_temp, color in zip(sorted(estimated_shift_factors.keys()), colors):
            ref_shift_factors = estimated_shift_factors[ref_temp]
            combined_log_freq = []
            combined_storage_modulus = []
            label_added = False
            for temp in sorted(temperatures):
                if temp in ref_shift_factors:
                    subset = data[data['Temperature'] == temp].copy()
                    subset['Adjusted Frequency (Hz)'] = subset['Frequency'] * ref_shift_factors[temp]
                    combined_log_freq.extend(np.log10(subset['Adjusted Frequency (Hz)']))
                    combined_storage_modulus.extend(subset['Storage Modulus'])
                    # ax.semilogx(subset['Adjusted Frequency (Hz)'], subset['Storage Modulus'], 'o', markersize=4, color=color)
                    if not label_added:
                        ax.semilogx(subset['Adjusted Frequency (Hz)'], subset['Storage Modulus'],
                            'o', markersize=4, color=color, label=f"{ref_temp}°C")
                        label_added = True
                    else:
                        ax.semilogx(subset['Adjusted Frequency (Hz)'], subset['Storage Modulus'],
                            'o', markersize=4, color=color)

            combined_log_freq = np.array(combined_log_freq)
            combined_storage_modulus = np.array(combined_storage_modulus)

            try:
                params, _ = curve_fit(storage_modulus_model,
                                      combined_log_freq, combined_storage_modulus,
                                      bounds=([1e-6, -100, -100, 1e-6], [a_fit, 100, 100, d_fit]),
                                      maxfev=9999999)
                # b_fit, c_fit = params[1], params[2]  
                a_fit_final, b_fit, c_fit, d_fit_final = params
                abcd_parameters.append((ref_temp, a_fit_final, b_fit, c_fit, d_fit_final))  
            except Exception as e:
                print(f"Curve fitting error for reference temp {ref_temp}: {e}")
                continue

            fitted_storage_modulus = storage_modulus_model(combined_log_freq, a_fit_final, b_fit, c_fit, d_fit_final)

            ax.plot(10**combined_log_freq, fitted_storage_modulus, color='black', linestyle='--', linewidth=1)

        ax.set_xlabel('Reduced Frequency (Hz)', fontname="Times New Roman", fontsize=20)
        ax.set_ylabel('Storage Modulus (MPa)', fontname="Times New Roman", fontsize=20)
        ax.set_title('Master Curve for All Reference Temperatures', fontname="Times New Roman", fontsize=20)
        ax.set_yscale('linear')
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontname("Times New Roman")
            label.set_fontsize(16)
        # ax.legend(fontsize=20, prop={"family": "Times New Roman"})

        ymin_val = np.min(combined_storage_modulus)
        ymax_val = np.max(combined_storage_modulus)
        if ymin_val < 0:
            lower_bound = 0
        else:
            lower_bound = ymin_val * 0
        ax.set_ylim(lower_bound, ymax_val * 1.1)

        add_watermark(ax, "NYU-ViscoMOD", on_axes=True)
        fig.tight_layout()

        fig.subplots_adjust(right=0.8)

        from matplotlib.lines import Line2D
        dummy_line = Line2D([], [], color='black', linestyle='--', linewidth=1, label='Fit Line')
        handles, labels = ax.get_legend_handles_labels()
        handles.insert(0, dummy_line)
        labels.insert(0, dummy_line.get_label())

        ax.legend(handles, labels,
            loc='center left',
            bbox_to_anchor=(1.01, 0.5),
            borderaxespad=0.,
            prop={"family": "Times New Roman", "size": 12})

        canvas_10.draw()  
        abcd_parameters_global = abcd_parameters  
        # print("ABCD Parameters:", abcd_parameters_global)  

    lbl_a = Label(top_frame, text="Upper Bound for 'a':", font=('Times New Roman', 12))
    lbl_a.pack(side="left", padx=5)
    scale_var_a = DoubleVar(value=a_fit_1)  
    a_slider = tk.Scale(top_frame, from_=0, to=2000, resolution=1, orient=HORIZONTAL, length=200,  
                        variable=scale_var_a, command=lambda val: update_graph(float(val), scale_var_d.get()))
    a_slider.pack(side="left", padx=5)

    btn_minus_a = Button(top_frame, text="-", font=('Times New Roman', 12), command=lambda: update_slider(scale_var_a, -1, lambda: update_graph(scale_var_a.get(), scale_var_d.get())))
    btn_minus_a.pack(side="left", padx=2)
    btn_plus_a = Button(top_frame, text="+", font=('Times New Roman', 12), command=lambda: update_slider(scale_var_a, 1, lambda: update_graph(scale_var_a.get(), scale_var_d.get())))
    btn_plus_a.pack(side="left", padx=2)

    lbl_d = Label(top_frame, text="Upper Bound for 'd':", font=('Times New Roman', 12))
    lbl_d.pack(side="left", padx=5)
    scale_var_d = DoubleVar(value=d_fit_1)  
    d_slider = tk.Scale(top_frame, from_=0, to=2000, resolution=1, orient=HORIZONTAL, length=200,  
                        variable=scale_var_d, command=lambda val: update_graph(scale_var_a.get(), float(val)))
    d_slider.pack(side="left", padx=5)

    btn_minus_d = Button(top_frame, text="-", font=('Times New Roman', 12), command=lambda: update_slider(scale_var_d, -1, lambda: update_graph(scale_var_a.get(), scale_var_d.get())))
    btn_minus_d.pack(side="left", padx=2)
    btn_plus_d = Button(top_frame, text="+", font=('Times New Roman', 12), command=lambda: update_slider(scale_var_d, 1, lambda: update_graph(scale_var_a.get(), scale_var_d.get())))
    btn_plus_d.pack(side="left", padx=2)

    canvas_10 = FigureCanvasTkAgg(fig, master=page_frame)
    canvas_10.draw()
    canvas_10.get_tk_widget().pack(fill="both", expand=True)

    update_graph(a_slider.get(), d_slider.get())


#############################################################
# Page 6: Modulus vs Strain Rates for Different Temperatures
#############################################################
# strain_rates_to_plot = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
strain_rates_to_plot = [1e-5, 1e-4, 1e-3, 1e-2]
varCheck = {}

fig_page11 = None
df_results_global = None

def elastic_modulus_vs_strain_rate_page(content_frame):
    global data, abcd_parameters_global, strain_rates_to_plot, varCheck, fig_page11

    if data is None or data.empty:
        messagebox.showerror("Error", "No data loaded!")
        return
    if not abcd_parameters_global:
        messagebox.showerror("Error", "No abcd_parameters available. Please complete fitting first.")
        return

    for widget in content_frame.winfo_children():
        widget.destroy()

    page_frame = Frame(content_frame)
    page_frame.pack(fill="both", expand=True)

    top_frame = Frame(page_frame)
    top_frame.pack(side="top", fill="x", padx=5, pady=5)

    right_frame = Frame(page_frame)
    right_frame.pack(side="right", fill="both", expand=True)

    def compute_data_results():
        results = []
        for params in abcd_parameters_global:
            ref_temp, a_val, b_val, c_val, d_val = params
            final_cumulative_integrals = []
            for rate in strain_rates_to_plot:
                time_min_rate = strain_min / rate
                time_max_rate = strain_max / rate
                time_range_rate = np.linspace(time_min_rate, time_max_rate, num_steps)
                E_t_time_range_rate = Etime_time_cycle(time_range_rate, 500, a_val, b_val, c_val, d_val)
                # if torsion_data is not None:
                #     Stress_history_rate = E_t_time_range_rate / 2.0 * rate
                # elif tension_data is not None:
                #     Stress_history_rate = E_t_time_range_rate * rate
                # else:
                #     messagebox.showerror("Error", "No valid data found for stress calculation.")
                #     return None
                if tension_data is not None:
                    Stress_history_rate = E_t_time_range_rate * rate
                else:
                    Stress_history_rate = E_t_time_range_rate / 2.0 * rate
                cumulative_integral = np.trapz(Stress_history_rate, time_range_rate) / strain_max
                final_cumulative_integrals.append(cumulative_integral)
            results.append([ref_temp] + final_cumulative_integrals)
        columns = ['Ref Temp (°C)'] + [f'Strain Rate {r} (1/s)' for r in strain_rates_to_plot]
        return pd.DataFrame(results, columns=columns)

    def update_plot():
        nonlocal right_frame
        global fig_page11, df_results_global

        for w in right_frame.winfo_children():
            w.destroy()
        
        if df_results_global is None:
            data_results = compute_data_results()
            if data_results is None:
                return
            df_results_global = data_results
        else:
            data_results = df_results_global

        fig = Figure(figsize=(12, 6))
        ax = fig.add_subplot(111)

        for i, row in data_results.iterrows():
            ref_temp = row['Ref Temp (°C)']
            # if temp_vars[ref_temp].get():
            y_values = row[1:].tolist()
            ax.plot(strain_rates_to_plot, y_values, marker='o', linestyle='-', label=f"{ref_temp}°C")
        ax.set_xscale('log')
        ax.set_xlabel('Strain Rate (1/s)', fontname="Times New Roman", fontsize=20)
        if tension_data is not None:
            ax.set_ylabel('Elastic Modulus (MPa)', fontname="Times New Roman", fontsize=20)
            ax.set_title('Elastic Modulus vs Strain Rate for Different Temperatures',
                         fontname="Times New Roman", fontsize=24)
        else:
            ax.set_ylabel('Shear Modulus (MPa)', fontname="Times New Roman", fontsize=20)
            ax.set_title('Shear Modulus vs Strain Rate for Different Temperatures',
                         fontname="Times New Roman", fontsize=24)

        from matplotlib.ticker import LogFormatterMathtext
        formatter = LogFormatterMathtext()
        ax.xaxis.set_major_formatter(formatter)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontname("Times New Roman")
            label.set_fontsize(16)
        
        y_data = data_results.iloc[:, 1:].to_numpy().flatten()
        ymin_val = np.min(y_data)
        ymax_val = np.max(y_data)
        if ymin_val < 0:
            lower_bound = 0
        else:
            lower_bound = ymin_val * 0
        ax.set_ylim(lower_bound, ymax_val * 1.1)
        ax.tick_params(axis='both', labelsize=16)
        # ax.legend(fontsize=15, prop={"family": "Times New Roman"})
        ax.grid(True)
        fig.tight_layout()

        fig.subplots_adjust(right=0.85)
        ax.legend(
            loc='center left',
            bbox_to_anchor=(1.01, 0.5),   
            borderaxespad=0.,
            prop={"family": "Times New Roman", "size": 12}
        )
        add_watermark(ax, "NYU-ViscoMOD", on_axes=True)
        fig_page11 = fig  

        canvas = FigureCanvasTkAgg(fig, master=right_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    btn_prev = Button(top_frame, text="Previous", font=('Times New Roman', 12),
                      command=lambda: show_page(master_curve_each_temp_page))
    btn_prev.pack(side="left", padx=5)

    btn_next = Button(top_frame, text="Next", font=('Times New Roman', 12),
                      command=lambda: show_page(modulus_vs_temperature_page))
    btn_next.pack(side="right", padx=5)
    btn_export = Button(top_frame, text="Export CSV", font=('Times New Roman', 12),
                        command=lambda: export_data_page11(compute_data_results()))
    btn_export.pack(side="right", padx=5)
    btn_save_graph = Button(top_frame, text="Save Graph", font=('Times New Roman', 12),
                            command=lambda: save_graph_page11(fig_page11))
    btn_save_graph.pack(side="right", padx=5)

    def export_data_page11(data_results):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All Files", "*.*")]
        )
        if file_path:
            data_results.to_csv(file_path, index=False)
            messagebox.showinfo("Info", f"Data exported to {file_path}")

    def save_graph_page11(fig):
        file_path = filedialog.asksaveasfilename(
            defaultextension=".jpeg",
            filetypes=[("JPEG files", "*.jpeg"), ("All Files", "*.*")]
        )
        if file_path and fig is not None:
            fig.savefig(file_path, dpi=300, format="jpeg", bbox_inches='tight')
            messagebox.showinfo("Info", f"Graph saved to {file_path}")

    update_plot()


############################################################
# Page 7: Modulus vs Temperature for Different Strain Rates
############################################################
df_results_global_page12 = None

def export_data_page12(data_results):
    file_path = filedialog.asksaveasfilename(
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All Files", "*.*")]
    )
    if file_path:
        data_results.to_csv(file_path, index=False)
        messagebox.showinfo("Info", f"Data exported to {file_path}")

def save_graph_page12(fig):
    file_path = filedialog.asksaveasfilename(
        defaultextension=".jpeg",
        filetypes=[("JPEG files", "*.jpeg"), ("All Files", "*.*")]
    )
    if file_path:
        fig.savefig(file_path, dpi=300, format="jpeg", bbox_inches='tight')
        messagebox.showinfo("Info", f"Graph saved to {file_path}")

def modulus_vs_temperature_page(content_frame):

    # clear content_frame
    for widget in content_frame.winfo_children():
        widget.destroy()
    global data, abcd_parameters_global
    if data is None or data.empty:
        messagebox.showerror("Error","No data loaded!")
        return
    if not abcd_parameters_global:
        messagebox.showerror("Error","No abcd_parameters available. Please complete Page 6 first.")
        return

    # Compute DataFrame 
    # results = []
    global df_results_global_page12
    if df_results_global_page12 is None:
        results = []
        for params in abcd_parameters_global:
            ref_temp, a_val, b_val, c_val, d_val = params
            final_cumulative_integrals = []
            for rate in strain_rates_to_plot:
                time_min_rate = strain_min / rate
                time_max_rate = strain_max / rate
                time_range_rate = np.linspace(time_min_rate, time_max_rate, num_steps)

                # Calculating E(t) using global Etime_time_cycle
                E_t_time_range_rate = Etime_time_cycle(time_range_rate, 500, a_val, b_val, c_val, d_val)
                if torsion_data is not None:
                    Stress_history_rate = E_t_time_range_rate/2 * rate
                elif tension_data is not None:
                    Stress_history_rate = E_t_time_range_rate * rate
                cumulative_integral = np.trapz(Stress_history_rate, time_range_rate) / strain_max
                final_cumulative_integrals.append(cumulative_integral)
            results.append([ref_temp] + final_cumulative_integrals)
        df_results_global_page12 = pd.DataFrame(
            results,
            columns=['Ref Temp (°C)'] + [f'Strain Rate {r} (1/s)' for r in strain_rates_to_plot]
        )
    data_results = df_results_global_page12
    # Top navigation
    top_frame = Frame(content_frame)
    top_frame.pack(side="top", fill="x", padx=5, pady=5)
    btn_prev = Button(
        top_frame, text="Previous", font=('Times New Roman', 12),
        command=lambda: show_page(elastic_modulus_vs_strain_rate_page)
    )
    btn_prev.pack(side="left", padx=5)
    btn_next = Button(top_frame, text="Next", font=('Times New Roman', 12),
                      command=lambda: show_page(elastic_modulus_3d_page))
    btn_next.pack(side="right", padx=5)

    btn_export = Button(top_frame, text="Export CSV", font=('Times New Roman', 12),
                        command=lambda: export_data_page12(data_results))
    btn_export.pack(side="right", padx=5)
    global fig_page12
    fig_page12 = None
    btn_save_graph = Button(top_frame, text="Save Graph", font=('Times New Roman', 12),
                            command=lambda: save_graph_page12(fig_page12))
    btn_save_graph.pack(side="right", padx=5)

    # Right plot area
    right_frame = Frame(content_frame)
    right_frame.pack(side="right", fill="both", expand=True)

    # Checkbox logic 
    varCheck = {}

    SUPERSCRIPT_MAP = {
        '0': '\u2070',
        '1': '\u00B9',
        '2': '\u00B2',
        '3': '\u00B3',
        '4': '\u2074',
        '5': '\u2075',
        '6': '\u2076',
        '7': '\u2077',
        '8': '\u2078',
        '9': '\u2079'
    }

    def to_superscript(exponent):

        sign_str = ''
        if exponent < 0:
            sign_str = '\u207B' 
            exponent = -exponent
        
        result = ''
        for ch in str(exponent):  
            result += SUPERSCRIPT_MAP[ch]
        return sign_str + result

    # Updating the plotting function 
    def update_modulus_vs_temp_plot():
        global fig_page12
        # Clear all child controls on the right (old canvas, etc.)
        for widget in right_frame.winfo_children():
            widget.destroy()

        # create Figure
        fig = Figure(figsize=(12, 6))
        ax = fig.add_subplot(111)

        # x axis
        x_vals = data_results['Ref Temp (°C)'].values

        # Iterate over all strain rates and plot if the check box is selected
        for rate in strain_rates_to_plot:
            #if varCheck[rate].get() == 1:
            col_name = f"Strain Rate {rate} (1/s)"
            y_vals = data_results[col_name].values
            # ax.plot(x_vals, y_vals, marker='s', linestyle='-', label=f"Strain Rate {rate} (1/s)")
            ax.plot(x_vals, y_vals, marker='s', linestyle='-', label=r"Strain Rate $10^{" + str(int(np.log10(rate))) + r"}$ (1/s)")

        ax.set_xlabel('Temperature (°C)', fontname="Times New Roman", fontsize=20)
        if tension_data is not None:
            ax.set_ylabel('Elastic Modulus (MPa)', fontname="Times New Roman", fontsize=20)
            ax.set_title('Elastic Modulus vs Temperature for Different Strain Rates',
                         fontname="Times New Roman", fontsize=24)
        elif torsion_data is not None:
            ax.set_ylabel('Shear Modulus (MPa)', fontname="Times New Roman", fontsize=20)
            ax.set_title('Shear Modulus vs Temperature for Different Strain Rates',
                         fontname="Times New Roman", fontsize=24)           
        # ax.set_ylim(bottom=0)

        # Adjust the scale font
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontname("Times New Roman")
            label.set_fontsize(16)
        
        y_values = data_results.iloc[:, 1:].values.flatten()
        ymin_val = np.min(y_values)
        ymax_val = np.max(y_values)
        if ymin_val < 0:
            lower_bound = 0
        else:
            lower_bound = ymin_val * 0
        ax.set_ylim(lower_bound, ymax_val * 1.1)
        ax.tick_params(axis='both', labelsize=16)
        # ax.legend(fontsize=15, prop={"family": "Times New Roman"})
        ax.grid(True)
        fig.tight_layout()

        fig.subplots_adjust(right=0.8)
        ax.legend(
            loc='center left',
            bbox_to_anchor=(1.01, 0.5),   
            borderaxespad=0.,
            prop={"family": "Times New Roman", "size": 12}
        )
        add_watermark(ax, "NYU-ViscoMOD", on_axes=True)
        fig_page12 = fig
        # Place the canvas inside right_frame
        canvas = FigureCanvasTkAgg(fig, master=right_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)

    update_modulus_vs_temp_plot()


###############################################
# Page 8: Elastic Modulus 3D Surface
###############################################
# use the variate to save figure
fig_page13 = None

def export_data_page13():
    global df_results_global
    if df_results_global is None or df_results_global.empty:
        messagebox.showinfo("Info", "No data to export!")
        return
    file_path = filedialog.asksaveasfilename(defaultextension=".csv",
                                             filetypes=[("CSV files", "*.csv"), ("All Files", "*.*")])
    if file_path:
        df_results_global.to_csv(file_path, index=False)
        messagebox.showinfo("Info", f"Data exported to {file_path}")

def save_graph_page13():
    global fig_page13
    if fig_page13 is None:
        messagebox.showerror("Error", "No figure to save!")
        return
    file_path = filedialog.asksaveasfilename(defaultextension=".jpeg",
                                             filetypes=[("JPEG files", "*.jpeg"), ("All Files", "*.*")])
    if file_path:
        fig_page13.savefig(file_path, dpi=300, format='jpeg', bbox_inches='tight')
        messagebox.showinfo("Info", f"Graph saved to {file_path}")

def elastic_modulus_3d_page(content_frame):

    for widget in content_frame.winfo_children():
        widget.destroy()
    
    global strain_rates, strain_rates_to_plot, log_strain_rates
    global df_results_global
    if df_results_global is None or df_results_global.empty:
        messagebox.showerror("Error", "No data loaded!")
        return

    page_frame = tk.Frame(content_frame)
    page_frame.pack(fill="both", expand=True)

    top_frame = tk.Frame(page_frame)
    top_frame.pack(side="top", fill="x", padx=5, pady=5)

    btn_prev = tk.Button(top_frame, text="Previous", font=('Times New Roman', 12),
                         command=lambda: show_page(modulus_vs_temperature_page))
    btn_prev.pack(side="left", padx=5)
    btn_export = Button(top_frame, text="Export CSV", font=('Times New Roman', 12),
                        command=export_data_page13)
    btn_export.pack(side="right", padx=5)
    btn_save_graph = Button(top_frame, text="Save Graph", font=('Times New Roman', 12),
                            command=save_graph_page13)
    btn_save_graph.pack(side="right", padx=5)

    # ============ 3D plot ============

    temperatures = df_results_global['Ref Temp (°C)'].values
    # strain_rates = np.array([0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
    strain_rates = strain_rates_to_plot
    elastic_modulus = df_results_global.iloc[:, 1:].values

    T_dense = np.linspace(temperatures.min(), temperatures.max(), 200)
    # S_dense = np.logspace(np.log10(strain_rates.min()), np.log10(strain_rates.max()), 200)
    S_dense = np.logspace(np.log10(min(strain_rates)), np.log10(max(strain_rates)), 200)
    T_mesh, S_mesh = np.meshgrid(T_dense, S_dense)

    original_points = np.array([(T, np.log10(S)) 
                                for T in temperatures 
                                for S in strain_rates])
    original_values = elastic_modulus.flatten()

    try:
        E_dense = griddata(original_points,
                           original_values,
                           (T_mesh, np.log10(S_mesh)),
                           method='cubic')
    except ValueError as e:
        if "different number of values and points" in str(e):
            messagebox.showerror("Error", "Different number of values and points. Please complete Page 5 first")
            return
        else:
            raise
    
    global fig_page13
    fig_page13 = Figure(figsize=(8,6))
    ax = fig_page13.add_subplot(111, projection='3d')

    surf = ax.plot_surface(T_mesh, 
                           np.log10(S_mesh), 
                           E_dense, 
                           cmap='rainbow')

    ax.set_xlabel('Temperature (°C)', fontname="Times New Roman", fontsize=18)
    ax.set_ylabel('Strain Rate (1/s)', fontname="Times New Roman", fontsize=18)
    if tension_data is not None:
        ax.set_zlabel('Elastic Modulus (MPa)', fontname="Times New Roman", fontsize=18)
    elif torsion_data is not None:
        ax.set_zlabel('Shear Modulus (MPa)', fontname="Times New Roman", fontsize=18)
    for label in ax.get_xticklabels() + ax.get_yticklabels() + ax.get_zticklabels():
        label.set_fontname("Times New Roman")
        label.set_fontsize(14)
    log_strain_rates = np.log10(strain_rates)
    ax.set_yticks(log_strain_rates)
    ax.set_yticklabels([f'$10^{{{int(rate)}}}$' for rate in log_strain_rates])
    # ax.set_yticklabels([f"10^{int(rate)}" for rate in log_strain_rates])
    for label in ax.get_yticklabels():
        label.set_fontname("Times New Roman")
        label.set_fontsize(14)
    
    cax = fig_page13.add_axes([0.86, 0.15, 0.03, 0.7])  # (left, bottom, width, height)
    color_bar = fig_page13.colorbar(surf, cax=cax, aspect=20)
    if tension_data is not None:
        color_bar.set_label('Elastic Modulus (MPa)', fontname="Times New Roman", fontsize=18)
    elif torsion_data is not None:
        color_bar.set_label('Shear Modulus (MPa)', fontname="Times New Roman", fontsize=18)
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray')
    for t in color_bar.ax.get_yticklabels():
        t.set_fontname('Times New Roman')
        t.set_fontsize(14)

    ax.view_init(elev=30, azim=340)

    # watermark
    add_watermark(ax, "NYU-ViscoMOD", on_axes=True)

    # put Figure in Tkinter 
    canvas = FigureCanvasTkAgg(fig_page13, master=page_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill="both", expand=True)


###############################################
# Sidebar navigation
###############################################
def setup_sidebar():
    global sidebar
    sidebar = Frame(root, width=150, bg="lightgray")
    sidebar.pack(side="left", fill="y")
    btn0 = Button(sidebar, text="Page 0\nLoad CSV", font=("Times New Roman", 12),
                  command=lambda: show_page(load_csv_page))
    btn0.pack(fill="x", padx=5, pady=6)
    btn1 = Button(sidebar, text="Page 1\nRaw Data", font=("Times New Roman", 12),
                  command=lambda: show_page(raw_data_page))
    btn1.pack(fill="x", padx=5, pady=6)
    btn2 = Button(sidebar, text="Page 2\n3D Surface", font=("Times New Roman", 12),
                  command=lambda: show_page(show_surface_plot))
    btn2.pack(fill="x", padx=5, pady=6)
    btn3 = Button(sidebar, text="Page 3\nMaster Curve", font=("Times New Roman", 12),
                  command=lambda: show_page(analysis_page))
    btn3.pack(fill="x", padx=5, pady=6)

    btn7 = Button(sidebar, text="Page 4\nCurve Fitting", font=("Times New Roman", 12),
                  command=go_to_curve_fitting_page)
    btn7.pack(fill="x", padx=5, pady=6)

    btn10 = Button(sidebar, text="Page 5\nMaster Curve Fit", font=("Times New Roman", 12),
                  command=lambda: show_page(master_curve_each_temp_page))
    btn10.pack(fill="x", padx=5, pady=6)
    btn11 = Button(sidebar, text="Page 6\nModulus vs Strain Rates", font=("Times New Roman", 12),
                  command=lambda: show_page(elastic_modulus_vs_strain_rate_page))
    btn11.pack(fill="x", padx=5, pady=6)
    btn12 = Button(sidebar, text="Page 7\nModulus vs Temperature", font=("Times New Roman", 12),
                  command=lambda: show_page(modulus_vs_temperature_page))
    btn12.pack(fill="x", padx=5, pady=6)
    btn13 = Button(sidebar, text="Page 8\nModulus 3D Surface", font=("Times New Roman", 12),
                  command=lambda: show_page(elastic_modulus_3d_page))
    btn13.pack(fill="x", padx=5, pady=6)

# def show_page(page_func):
#     global content_frame
#     if page_func != load_csv_page:
#         if not hasattr(root, 'sidebar_created'):
#             setup_sidebar()
#             root.sidebar_created = True
#     if 'content_frame' not in globals():
#         content_frame = Frame(root)
#         content_frame.pack(side="right", fill="both", expand=True)
#     for widget in content_frame.winfo_children():
#         widget.destroy()
#     page_func(content_frame)
# def show_page(page_func):
#     global content_frame
    
#     for widget in content_frame.winfo_children():
#         widget.destroy()
#     page_func(content_frame)

def show_page(page_func):
    global content_frame, sidebar
    if page_func == load_csv_page:
        sidebar.pack_forget()
    else:
        sidebar.pack(side="left", fill="y")
    for widget in content_frame.winfo_children():
        widget.destroy()
    page_func(content_frame)

###############################################
# Entrance
###############################################
def main():
    # Creating a content area
    global content_frame
    setup_sidebar()
    content_frame = Frame(root)
    content_frame.pack(side="right", fill="both", expand=True)
    # By default, page 1 is displayed
    show_page(load_csv_page)
    root.mainloop()

if __name__ == "__main__":
    main()