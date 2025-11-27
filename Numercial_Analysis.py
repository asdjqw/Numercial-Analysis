import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib
matplotlib.use('TkAgg')

plt.rcParams["font.family"] = ["SimHei"]
plt.rcParams["axes.unicode_minus"] = False 
plt.rcParams['font.size'] = 12

def draw_function():
    root = tk.Tk()
    root.title("数值分析实验——迭代法解方程")
    root.geometry("1500x800")
    func_var = tk.StringVar(value="x^2 - 2")
    x_start_var = tk.StringVar(value="-5")
    x_end_var = tk.StringVar(value="5")
    x0_var = tk.StringVar(value="1")
    x1_var = tk.StringVar(value="2")
    iter_var = tk.StringVar(value="20")
    iter_func_var = tk.StringVar(value="x - (x^2 - 2)/(2*x)")
    
    main_frame = ttk.Frame(root, padding="30")
    main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
    
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    main_frame.columnconfigure(0, weight=1)
    main_frame.columnconfigure(1, weight=3)
    main_frame.columnconfigure(2, weight=1)
    main_frame.rowconfigure(6, weight=1)
    main_frame.rowconfigure(8, weight=2)
    
    #函数表达式输入
    ttk.Label(main_frame, text="原函数表达式（自变量x）：", font=("Arial", 14)).grid(row=0, column=0, sticky=tk.W, pady=8)
    func_entry = ttk.Entry(main_frame, textvariable=func_var, width=40, font=("Arial", 14))
    func_entry.grid(row=0, column=1, sticky=(tk.W, tk.E), pady=8, padx=5)
    ttk.Label(main_frame, text="示例：x^2+sin(2*x) 支持sin、cos、exp、log、sqrt", font=("Arial", 11), foreground="gray").grid(row=0, column=2, sticky=tk.W, pady=8)
    
    #绘图区间输入
    ttk.Label(main_frame, text="绘图x起始值：", font=("Arial", 14)).grid(row=1, column=0, sticky=tk.W, pady=8)
    x_start_entry = ttk.Entry(main_frame, textvariable=x_start_var, width=20, font=("Arial", 14))
    x_start_entry.grid(row=1, column=1, sticky=tk.W, pady=8, padx=5)
    ttk.Label(main_frame, text="绘图x结束值：", font=("Arial", 14)).grid(row=1, column=1, sticky=tk.E, pady=8, padx=20)
    x_end_entry = ttk.Entry(main_frame, textvariable=x_end_var, width=20, font=("Arial", 14))
    x_end_entry.grid(row=1, column=2, sticky=tk.W, pady=8, padx=5)
    
    #迭代参数输入
    ttk.Label(main_frame, text="迭代初值x0：", font=("Arial", 14)).grid(row=2, column=0, sticky=tk.W, pady=8)
    x0_entry = ttk.Entry(main_frame, textvariable=x0_var, width=20, font=("Arial", 14))
    x0_entry.grid(row=2, column=1, sticky=tk.W, pady=8, padx=5)
    ttk.Label(main_frame, text="迭代初值x1（仅弦截法）：", font=("Arial", 14)).grid(row=2, column=1, sticky=tk.E, pady=8, padx=20)
    x1_entry = ttk.Entry(main_frame, textvariable=x1_var, width=20, font=("Arial", 14))
    x1_entry.grid(row=2, column=2, sticky=tk.W, pady=8, padx=5)
    ttk.Label(main_frame, text="迭代次数：", font=("Arial", 14)).grid(row=3, column=0, sticky=tk.W, pady=8)
    iter_entry = ttk.Entry(main_frame, textvariable=iter_var, width=20, font=("Arial", 14))
    iter_entry.grid(row=3, column=1, sticky=tk.W, pady=8, padx=5)
    
    #迭代函数输入（艾特肯法、普通迭代法）
    ttk.Label(main_frame, text="迭代函数（艾特肯/普通迭代）：", font=("Arial", 14)).grid(row=4, column=0, sticky=tk.W, pady=8)
    iter_func_entry = ttk.Entry(main_frame, textvariable=iter_func_var, width=40, font=("Arial", 14))
    iter_func_entry.grid(row=4, column=1, sticky=(tk.W, tk.E), pady=8, padx=5)
    ttk.Label(main_frame, text="示例：x - (x^2-2)/(2*x)（牛顿迭代等价）", font=("Arial", 11), foreground="gray").grid(row=4, column=2, sticky=tk.W, pady=8)
    buttonf = ttk.Frame(main_frame)
    buttonf.grid(row=5, column=0, columnspan=3, pady=20)
    style = ttk.Style()
    style.configure("TButton", font=("Arial", 14), padding=10)
    style.configure("Accent.TButton", font=("Arial", 14, "bold"), padding=10, background="#4CAF50")

    #调整按钮间距
    ttk.Button(buttonf, text="1. 采样绘图+根区间检测", command=lambda: run_option(1), style="Accent.TButton").grid(row=0, column=0, padx=10, pady=8)
    ttk.Button(buttonf, text="2. 牛顿迭代法求根", command=lambda: run_option(2), style="Accent.TButton").grid(row=0, column=1, padx=10, pady=8)
    ttk.Button(buttonf, text="3. 双点弦截法求根", command=lambda: run_option(3), style="Accent.TButton").grid(row=0, column=2, padx=10, pady=8)
    ttk.Button(buttonf, text="4. 艾特肯法求根", command=lambda: run_option(4), style="Accent.TButton").grid(row=0, column=3, padx=10, pady=8)
    ttk.Button(buttonf, text="5. 普通迭代法求根", command=lambda: run_option(5), style="Accent.TButton").grid(row=0, column=4, padx=10, pady=8)
    ttk.Button(buttonf, text="6. 退出程序", command=root.quit, style="Accent.TButton").grid(row=0, column=5, padx=10, pady=8)
    ttk.Label(main_frame, text="迭代结果输出：", font=("Arial", 16, "bold")).grid(row=6, column=0, sticky=tk.W, pady=10)
    result_text = scrolledtext.ScrolledText(main_frame, width=120, height=25, font=("Consolas", 12))
    result_text.grid(row=7, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=10)
    fig, ax = plt.subplots(figsize=(12, 7))
    canvas = FigureCanvasTkAgg(fig, master=main_frame)
    canvas.draw()
    toolbar_frame = ttk.Frame(main_frame)
    toolbar_frame.grid(row=9, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
    toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
    toolbar.pack(side=tk.TOP, fill=tk.X)
    toolbar.update()
    for widget in toolbar.winfo_children():
        if isinstance(widget, tk.Button):
            widget.config(font=("Arial", 10))
    canvas.get_tk_widget().grid(row=8, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=10)
    map_1 = {
        'sin': 'np.sin',
        'cos': 'np.cos',
        'tan': 'np.tan',
        'exp': 'np.exp',
        'log10': 'np.log10',
        'log': 'np.log',
        'sqrt': 'np.sqrt',
        'pi': 'np.pi',
        '^': '**'
    }
    map_2 = {
        'sin': 'sp.sin',
        'cos': 'sp.cos',
        'tan': 'sp.tan',
        'exp': 'sp.exp',
        'log10': 'sp.log10',
        'log': 'sp.log',
        'sqrt': 'sp.sqrt',
        'pi': 'sp.pi',
        '^': '**'
    }
    
    def print_result(msg, clear=False):
        """在滚动文本框中显示结果，支持清空和换行"""
        if clear:
            result_text.delete(1.0, tk.END)
        result_text.insert(tk.END, msg + "\n")
        result_text.see(tk.END)
    
    def dofancf(func_str, x_value):
        """计算原函数值，避免x未定义错误"""
        func_processed = func_str
        for old, new in map_1.items():
            func_processed = func_processed.replace(old, new)
        
        try:
            result = eval(func_processed, {'np': np, 'x': float(x_value)})
            return result
        except Exception as e:
            try:
                func = eval(f"lambda x: {func_processed}", {'np': np})
                return func(float(x_value))
            except Exception as e2:
                print_result(f"计算原函数值失败：{str(e2)}")
                return None
    def to_sympy_num(value):
        """将值转换为sympy数值类型，避免int/float调用evalf()报错"""
        if isinstance(value, (int, float)):
            return sp.S(value)
        return value
    
    #执行
    def run_option(flag):
        f0 = func_var.get().strip()
        if not f0:
            messagebox.showwarning("警告", "请输入函数表达式！")
            return
        print_result(f"=== 执行操作：{['', '采样绘图+根区间检测', '牛顿迭代法', '双点弦截法', '艾特肯法', '普通迭代法'][flag]} ===", clear=True)
        print_result(f"使用函数：y = {f0}")
        
        try:
            if flag == 1:
                f0_original = f0
                f_processed = f0
                for old, new in map_1.items():
                    f_processed = f_processed.replace(old, new)
                try:
                    x_start = float(x_start_var.get().strip())
                    x_end = float(x_end_var.get().strip())
                    if x_start >= x_end:
                        print_result("错误：起始值必须小于结束值")
                        messagebox.showerror("错误", "起始值必须小于结束值")
                        return
                except ValueError:
                    print_result("错误：x区间必须输入数字")
                    messagebox.showerror("错误", "x区间必须输入数字")
                    return
                
                N = 10000
                x_vals = np.linspace(x_start, x_end, N)
                try:
                    y_vals = eval(f_processed, {'np': np, 'x': x_vals})
                except Exception as e:
                    print_result(f"函数表达式错误：{e}")
                    messagebox.showerror("错误", f"函数表达式错误：{e}")
                    return
                
                roots = []
                sign_changes = []
                for i in range(N - 1):
                    y1, y2 = y_vals[i], y_vals[i + 1]
                    if np.isnan(y1) or np.isnan(y2):
                        continue
                    if y1 == 0:
                        sign_changes.append(x_vals[i])
                    elif y1 * y2 < 0:
                        roots.append((x_vals[i], x_vals[i + 1]))
                        sign_changes.append((x_vals[i] + x_vals[i + 1]) / 2)
                if roots:
                    print_result(f"\n检测到 {len(roots)} 个可能包含实根的区间：")
                    print_result("-" * 60)
                    print_result(f"{'序号':<8}{'区间':<35}{'疑似根位置':<25}")
                    print_result("-" * 60)
                    for idx, (a, b) in enumerate(roots, 1):
                        mid = (a + b) / 2
                        print_result(f"{idx:<8}[{a:.6f}, {b:.6f}]{mid:.10f}")
                else:
                    print_result("\n未检测到符号变化，可能无实根或根为重根（需进一步分析）")
                ax.clear()
                ax.plot(x_vals, y_vals, color='blue', linewidth=2.0, label=f"y = {f0_original}")
                
                if sign_changes:
                    ax.scatter(sign_changes, [0]*len(sign_changes), color='red', s=50, zorder=5,
                               label='疑似根位置', marker='o')
                
                ax.axhline(y=0, color='black', linewidth=1.0, alpha=0.8)
                ax.axvline(x=0, color='black', linewidth=1.0, alpha=0.8)
                ax.set_xlabel('x', fontsize=14, fontweight='bold')
                ax.set_ylabel('y', fontsize=14, fontweight='bold')
                ax.set_title('函数图像与自动根区间检测', fontsize=16, fontweight='bold', pad=20)
                ax.grid(True, linestyle='--', alpha=0.8, linewidth=1.0)
                ax.legend(fontsize=12)
                plt.tight_layout()
                canvas.draw()
            
            elif flag == 2:
                f_processed = f0
                for old, new in map_2.items():
                    f_processed = f_processed.replace(old, new)
                
                x = sp.Symbol('x')
                try:
                    f0_sympy = eval(f_processed)
                except Exception as e:
                    print_result(f"函数表达式错误：{e}")
                    messagebox.showerror("错误", f"函数表达式错误：{e}")
                    return
                try:
                    f1 = sp.diff(f0_sympy, x)
                except Exception as e:
                    print_result(f"求导失败：{e}")
                    messagebox.showerror("错误", f"求导失败：{e}")
                    return
                x0_str = x0_var.get().strip()
                x0_processed = x0_str
                for old, new in map_2.items():
                    x0_processed = x0_processed.replace(old, new)
                try:
                    x0 = eval(x0_processed, {'sp': sp, 'np': np})
                    x0 = to_sympy_num(x0)
                except Exception as e:
                    print_result(f"初值输入错误：{e}")
                    messagebox.showerror("错误", f"初值输入错误：{e}")
                    return
                
                try:
                    time = int(iter_var.get().strip())
                    if time <= 0:
                        print_result("错误：迭代次数必须为正整数")
                        messagebox.showerror("错误", "迭代次数必须为正整数")
                        return
                except ValueError:
                    print_result("错误：迭代次数必须为整数")
                    messagebox.showerror("错误", "迭代次数必须为整数")
                    return
                print_result(f"\n牛顿迭代法参数：")
                print_result(f"初值x0 = {x0.evalf(10)}")
                print_result(f"迭代次数 = {time}")
                print_result(f"收敛精度 = 1e-10")
                print_result("\n迭代过程：")
                print_result("-" * 90)
                print_result(f"{'迭代次数':<12}{'当前值x':<28}{'函数值f(x)':<28}{'导数f’(x)':<22}")
                print_result("-" * 90)
                
                converged = False
                for i in range(time):
                    f1_val = f1.subs(x, x0).evalf()
                    f0_val = f0_sympy.subs(x, x0).evalf()
                    print_result(f"{i+1:<12}{x0.evalf(10):<28}{f0_val.evalf(10):<28}{f1_val.evalf(10):<22}")
                    
                    if abs(f1_val) < 1e-10:
                        print_result("\n警告：导数接近0，无法继续迭代")
                        messagebox.showwarning("警告", "导数接近0，无法继续迭代")
                        break
                    
                    x1 = x0 - f0_val / f1_val
                    x1 = to_sympy_num(x1)
                    if abs(x1.evalf() - x0.evalf()) < 1e-10:
                        print_result("-" * 90)
                        print_result(f"\n迭代{i+1}次收敛！")
                        print_result(f"最终根值：{x1.evalf(10)}")
                        print_result(f"最终函数值：{f0_sympy.subs(x, x1).evalf(10)}")
                        converged = True
                        break 
                    x0 = x1
                    if i >= 5000:
                        print_result("\n警告：超过最大迭代次数（5000），请重新选择初值")
                        messagebox.showwarning("警告", "超过最大迭代次数，请重新选择初值")
                        break
                if not converged and time <= 5000:
                    print_result("-" * 90)
                    print_result(f"\n完成{time}次迭代（未完全收敛）")
                    print_result(f"最终迭代值：{x0.evalf(10)}")
                    print_result(f"最终函数值：{f0_sympy.subs(x, x0).evalf(10)}")
            
            elif flag == 3:
                f_processed = f0
                for old, new in map_2.items():
                    f_processed = f_processed.replace(old, new)
                x = sp.Symbol('x')
                try:
                    f0_sympy = eval(f_processed)
                except Exception as e:
                    print_result(f"函数表达式错误：{e}")
                    messagebox.showerror("错误", f"函数表达式错误：{e}")
                    return
                x0_str = x0_var.get().strip()
                x1_str = x1_var.get().strip()
                x0_processed = x0_str
                x1_processed = x1_str
                for old, new in map_2.items():
                    x0_processed = x0_processed.replace(old, new)
                    x1_processed = x1_processed.replace(old, new)
                try:
                    x0 = eval(x0_processed, {'sp': sp, 'np': np})
                    x1 = eval(x1_processed, {'sp': sp, 'np': np})
                    x0 = to_sympy_num(x0)
                    x1 = to_sympy_num(x1)
                except Exception as e:
                    print_result(f"初值输入错误：{e}")
                    messagebox.showerror("错误", f"初值输入错误：{e}")
                    return
                try:
                    time = int(iter_var.get().strip())
                    if time <= 0:
                        print_result("错误：迭代次数必须为正整数")
                        messagebox.showerror("错误", "迭代次数必须为正整数")
                        return
                except ValueError:
                    print_result("错误：迭代次数必须为整数")
                    messagebox.showerror("错误", "迭代次数必须为整数")
                    return
                print_result(f"\n双点弦截法参数：")
                print_result(f"初值x0 = {x0.evalf(10)}, x1 = {x1.evalf(10)}")
                print_result(f"迭代次数 = {time}")
                print_result(f"收敛精度 = 1e-10")
                print_result("\n迭代过程：")
                print_result("-" * 80)
                print_result(f"{'迭代次数':<12}{'x_prev':<25}{'x_curr':<25}{'x_new':<25}")
                print_result("-" * 80)
                converged = False
                for i in range(time):
                    f0_val = f0_sympy.subs(x,x0).evalf()
                    f1_val = f0_sympy.subs(x,x1).evalf()
                    
                    if abs(f1_val-f0_val) < 1e-10:
                        print_result("\n警告：分母接近0，无法继续迭代")
                        messagebox.showwarning("警告", "分母接近0，无法继续迭代")
                        break
                    x2 = x1 - f1_val * (x1 - x0) / (f1_val - f0_val)
                    x2 = to_sympy_num(x2)
                    print_result(f"{i+1:<12}{x0.evalf(10):<25}{x1.evalf(10):<25}{x2.evalf(10):<25}")
                    if abs(x2.evalf() - x1.evalf()) < 1e-10:
                        print_result("-" * 80)
                        print_result(f"\n迭代{i+1}次收敛！")
                        print_result(f"最终根值：{x2.evalf(10)}")
                        print_result(f"最终函数值：{f0_sympy.subs(x, x2).evalf(10)}")
                        converged = True
                        break
                    x0 = x1
                    x1 = x2
                    if i >= 5000:
                        print_result("\n警告：超过最大迭代次数（5000），请重新选择初值")
                        messagebox.showwarning("警告", "超过最大迭代次数，请重新选择初值")
                        break
                if not converged and time <= 5000:
                    print_result("-" * 80)
                    print_result(f"\n完成{time}次迭代（未完全收敛）")
                    print_result(f"最终迭代值：{x2.evalf(10)}")
                    print_result(f"最终函数值：{f0_sympy.subs(x, x2).evalf(10)}")
            
            elif flag == 4:
                fp = iter_func_var.get().strip()
                if not fp:
                    print_result("错误：请输入迭代函数")
                    messagebox.showerror("错误", "请输入迭代函数")
                    return
                print_result(f"迭代函数：φ(x) = {fp}")
                fp_processed = fp
                for old, new in map_2.items():
                    fp_processed = fp_processed.replace(old, new)
                x = sp.Symbol('x')
                try:
                    f0_sympy = eval(fp_processed)
                except Exception as e:
                    print_result(f"迭代函数表达式错误：{e}")
                    messagebox.showerror("错误", f"迭代函数表达式错误：{e}")
                    return
                x0_str = x0_var.get().strip()
                x0_processed = x0_str
                for old, new in map_2.items():
                    x0_processed = x0_processed.replace(old, new)
                try:
                    x0 = eval(x0_processed, {'sp': sp, 'np': np})
                    x0 = to_sympy_num(x0)
                except Exception as e:
                    print_result(f"初值输入错误：{e}")
                    messagebox.showerror("错误", f"初值输入错误：{e}")
                    return
                try:
                    time = int(iter_var.get().strip())
                    if time <= 0:
                        print_result("错误：迭代次数必须为正整数")
                        messagebox.showerror("错误", "迭代次数必须为正整数")
                        return
                except ValueError:
                    print_result("错误：迭代次数必须为整数")
                    messagebox.showerror("错误", "迭代次数必须为整数")
                    return
                print_result(f"\n艾特肯法参数：")
                print_result(f"初值x0 = {x0.evalf(10)}")
                print_result(f"迭代次数 = {time}")
                print_result(f"收敛精度 = 1e-10")
                print_result("\n迭代过程：")
                print_result("-" * 100)
                print_result(f"{'迭代次数':<12}{'x0':<25}{'y1=φ(x0)':<25}{'z1=φ(y1)':<25}{'x1(艾特肯)':<25}")
                print_result("-" * 100)
                
                converged = False
                for i in range (time):
                    y1 = f0_sympy.subs(x,x0).evalf()
                    z1 = f0_sympy.subs(x,y1).evalf()
                    y1 = to_sympy_num(y1)
                    z1 = to_sympy_num(z1)
                    if abs(y1.evalf() - x0.evalf()) < 1e-10:
                        print_result("\n警告：分母接近0，无法继续迭代")
                        messagebox.showwarning("警告", "分母接近0，无法继续迭代")
                        break
                    
                    g = (z1.evalf() - y1.evalf()) / (y1.evalf() - x0.evalf())
                    x1 = (1 / (1 - g)) * (y1.evalf() - g * x0.evalf())
                    x1 = to_sympy_num(x1)
                    print_result(f"{i+1:<12}{x0.evalf(10):<25}{y1.evalf(10):<25}{z1.evalf(10):<25}{x1.evalf(10):<25}")
                    
                    if abs(x1.evalf() - x0.evalf()) < 1e-10:
                        print_result("-" * 100)
                        print_result(f"\n迭代{i+1}次收敛！")
                        print_result(f"最终根值：{x1.evalf(10)}")
                        orig_func_val = dofancf(f0, x1.evalf())
                        if orig_func_val is not None:
                            print_result(f"原函数值：{orig_func_val:.10f}")
                        converged = True
                        break
                    x0 = x1
                    if i >= 5000:
                        print_result("\n警告：超过最大迭代次数（5000），请重新选择初值")
                        messagebox.showwarning("警告", "超过最大迭代次数，请重新选择初值")
                        break
                if not converged and time <= 5000:
                    print_result("-" * 100)
                    print_result(f"\n完成{time}次迭代（未完全收敛）")
                    print_result(f"最终迭代值：{x1.evalf(10)}")
                    orig_func_val = dofancf(f0, x1.evalf())
                    if orig_func_val is not None:
                        print_result(f"原函数值：{orig_func_val:.10f}")
            
            elif flag == 5:
                fp = iter_func_var.get().strip()
                if not fp:
                    print_result("错误：请输入迭代函数")
                    messagebox.showerror("错误", "请输入迭代函数")
                    return
                
                print_result(f"迭代函数：φ(x) = {fp}")
                fp_processed = fp
                for old, new in map_2.items():
                    fp_processed = fp_processed.replace(old, new)
                x = sp.Symbol('x')
                try:
                    f0_sympy = eval(fp_processed)
                except Exception as e:
                    print_result(f"迭代函数表达式错误：{e}")
                    messagebox.showerror("错误", f"迭代函数表达式错误：{e}")
                    return
                x0_str = x0_var.get().strip()
                x0_processed = x0_str
                for old, new in map_2.items():
                    x0_processed = x0_processed.replace(old, new)
                try:
                    x0 = eval(x0_processed, {'sp': sp, 'np': np})
                    x0 = to_sympy_num(x0)
                except Exception as e:
                    print_result(f"初值输入错误：{e}")
                    messagebox.showerror("错误", f"初值输入错误：{e}")
                    return
                try:
                    time = int(iter_var.get().strip())
                    if time <= 0:
                        print_result("错误：迭代次数必须为正整数")
                        messagebox.showerror("错误", "迭代次数必须为正整数")
                        return
                except ValueError:
                    print_result("错误：迭代次数必须为整数")
                    messagebox.showerror("错误", "迭代次数必须为整数")
                    return
                print_result(f"\n普通迭代法参数：")
                print_result(f"初值x0 = {x0.evalf(10)}")
                print_result(f"迭代次数 = {time}")
                print_result(f"收敛精度 = 1e-10")
                print_result("\n迭代过程：")
                print_result("-" * 70)
                print_result(f"{'迭代次数':<12}{'当前值x':<28}{'下一次值φ(x)':<28}")
                print_result("-" * 70)
                converged = False
                for i in range (time):
                    x1 = f0_sympy.subs(x,x0).evalf()
                    x1 = to_sympy_num(x1)
                    print_result(f"{i+1:<12}{x0.evalf(10):<28}{x1.evalf(10):<28}")
                    
                    if abs(x1.evalf() - x0.evalf()) < 1e-10:
                        print_result("-" * 70)
                        print_result(f"\n迭代{i+1}次收敛！")
                        print_result(f"最终根值：{x1.evalf(10)}")

                        orig_func_val = dofancf(f0, x1.evalf())
                        if orig_func_val is not None:
                            print_result(f"原函数值：{orig_func_val:.10f}")
                        converged = True
                        break
                    
                    x0 = x1
                    if i >= 5000:
                        print_result("\n警告：超过最大迭代次数（5000），请重新选择初值")
                        messagebox.showwarning("警告", "超过最大迭代次数，请重新选择初值")
                        break
                
                if not converged and time <= 5000:
                    print_result("-" * 70)
                    print_result(f"\n完成{time}次迭代（未完全收敛）")
                    print_result(f"最终迭代值：{x1.evalf(10)}")
                    orig_func_val = dofancf(f0, x1.evalf())
                    if orig_func_val is not None:
                        print_result(f"原函数值：{orig_func_val:.10f}")
        
        except Exception as e:
            print_result(f"程序异常：{str(e)}")
            messagebox.showerror("异常", f"程序异常：{str(e)}")
    
    root.mainloop()

if __name__ == "__main__":
    draw_function()
