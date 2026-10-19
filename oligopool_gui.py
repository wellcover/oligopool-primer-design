#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oligopool_gui.py - Oligopool 引物设计工具 (图形界面)

工作流分页: ①设计引物 → ②结果(表格/复制) → ③校验引物 → ④生成合成序列 → 日志
依赖同目录的 oligopool_primer_design.py 与 build_oligos.py, 结果与命令行一致。
"""
import contextlib
import os
import queue
import threading
import time
import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox
import webbrowser

import oligopool_primer_design as opd
import build_oligos as bo

APP_TITLE = "Oligopool 引物设计工具"
APP_VERSION = "1.1"
REPO_URL = "https://github.com/wellcover/oligopool-primer-design"
FONT_UI = ("Microsoft YaHei UI", 10)
FONT_UI_SMALL = ("Microsoft YaHei UI", 9)
FONT_MONO = ("Consolas", 10)
C_PASS = "#1a7f37"
C_FAIL = "#c62828"
C_MUTED = "#666666"


class Tooltip(object):
    # 悬浮提示
    def __init__(self, widget, text, delay=500):
        self.widget, self.text, self.delay = widget, text, delay
        self.tip, self.after_id = None, None
        widget.bind("<Enter>", self._schedule)
        widget.bind("<Leave>", self._schedule_hide)

    def _schedule(self, _e):
        self.after_id = self.widget.after(self.delay, self._show)

    def _schedule_hide(self, _e):
        if self.after_id:
            self.widget.after_cancel(self.after_id)
        self._hide()

    def _show(self):
        if self.tip:
            return
        x = self.widget.winfo_rootx() + 12
        y = self.widget.winfo_rooty() + self.widget.winfo_height() + 6
        self.tip = tw = tk.Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry("+%d+%d" % (x, y))
        tk.Label(tw, text=self.text, justify="left", bg="#ffffe0",
                 fg="#222222", font=FONT_UI_SMALL, relief="solid", bd=1,
                 padx=6, pady=3).pack()

    def _hide(self):
        if self.tip:
            self.tip.destroy()
            self.tip = None


class QueueWriter(object):
    # 把 print 输出转发到 GUI 队列
    def __init__(self, q):
        self.q = q

    def write(self, s):
        if s:
            self.q.put(s)

    def flush(self):
        pass


class App(object):
    def __init__(self, root):
        self.root = root
        root.title("%s v%s" % (APP_TITLE, APP_VERSION))
        root.minsize(980, 720)
        root.geometry("1080x780")
        self.q = queue.Queue()
        self.job = None
        self.job_done = None
        self.t_start = None
        self.busy = False
        self._last_outputs = []
        self._check_rows = []
        self._build_style()
        self._build_menu()
        self._build_widgets()
        self._set_busy(False, "就绪")
        self.root.after(120, self._poll)

    def _build_style(self):
        self.root.option_add("*Font", FONT_UI)
        style = ttk.Style(self.root)
        for theme in ("vista", "clam"):
            try:
                style.theme_use(theme)
                break
            except tk.TclError:
                continue
        style.configure("Primary.TButton", font=("Microsoft YaHei UI", 10, "bold"))
        style.configure("Step.TLabelframe.Label", font=("Microsoft YaHei UI", 10, "bold"))
        style.configure("Hint.TLabel", foreground=C_MUTED, font=FONT_UI_SMALL)
        style.configure("Status.TLabel", font=FONT_UI_SMALL)

    def _build_menu(self):
        m = tk.Menu(self.root)
        fm = tk.Menu(m, tearoff=0)
        fm.add_command(label="打开输出文件夹", command=self.open_outdir)
        fm.add_separator()
        fm.add_command(label="退出", command=self.root.destroy)
        m.add_cascade(label="文件", menu=fm)
        hm = tk.Menu(m, tearoff=0)
        hm.add_command(label="使用说明", command=self.show_help)
        hm.add_command(label="项目主页 (GitHub)", command=lambda: webbrowser.open(REPO_URL))
        hm.add_command(label="关于", command=self.show_about)
        m.add_cascade(label="帮助", menu=hm)
        self.root.config(menu=m)

    def _build_widgets(self):
        banner = ttk.Label(self.root, style="Hint.TLabel",
                           text="工作流: ① 设计引物 → ② 在结果中挑选两条 → ③ 可校验安全性 → ④ 填入 F/R 生成送合成序列      (k 阈值自动选择; 全部检查含反向互补形式)")
        banner.pack(fill="x", padx=10, pady=(6, 0))

        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill="both", expand=True, padx=10, pady=8)

        self.tab_design = ttk.Frame(self.nb)
        self.tab_results = ttk.Frame(self.nb)
        self.tab_check = ttk.Frame(self.nb)
        self.tab_build = ttk.Frame(self.nb)
        self.tab_log = ttk.Frame(self.nb)
        self.nb.add(self.tab_design, text=" ① 设计引物 ")
        self.nb.add(self.tab_results, text=" ② 结果 ")
        self.nb.add(self.tab_check, text=" ③ 校验引物 ")
        self.nb.add(self.tab_build, text=" ④ 生成合成序列 ")
        self.nb.add(self.tab_log, text=" 日志 ")

        self._build_design_tab()
        self._build_results_tab()
        self._build_check_tab()
        self._build_build_tab()
        self._build_log_tab()

        sb = ttk.Frame(self.root, style="Status.TLabel")
        sb.pack(fill="x", side="bottom", padx=10, pady=2)
        self.progress = ttk.Progressbar(sb, mode="indeterminate", length=140)
        self.progress.pack(side="left", padx=(0, 8))
        self.status = ttk.Label(sb, text="就绪", style="Status.TLabel", width=60, anchor="w")
        self.status.pack(side="left")
        self.elapsed = ttk.Label(sb, text="", style="Status.TLabel", anchor="e")
        self.elapsed.pack(side="right")


    # ---- 页1: 设计 ----
    def _build_design_tab(self):
        pane = ttk.PanedWindow(self.tab_design, orient="horizontal")
        pane.pack(fill="both", expand=True, padx=6, pady=6)

        left = ttk.Labelframe(pane, text=" ① 输入 FASTA 文件 ", style="Step.TLabelframe")
        pane.add(left, weight=1)
        self.file_list = tk.Listbox(left, height=8, selectmode="extended",
                                    font=FONT_MONO, activestyle="dotbox")
        self.file_list.pack(fill="both", expand=True, padx=(8, 2), pady=6)
        fb = ttk.Frame(left)
        fb.pack(fill="x", padx=8, pady=(0, 6))
        ttk.Button(fb, text="添加文件…", command=self.add_files).pack(side="left", padx=2)
        ttk.Button(fb, text="移除所选", command=self.remove_sel).pack(side="left", padx=2)
        ttk.Button(fb, text="清空", command=lambda: self.file_list.delete(0, "end")).pack(side="left", padx=2)

        mf = ttk.Labelframe(left, text=" 多文件模式 ", style="Step.TLabelframe")
        mf.pack(fill="x", padx=8, pady=(0, 8))
        self.mode = tk.StringVar(value="separate")
        ttk.Radiobutton(mf, text="模式1 · 每个文件单独一张芯片(依次独立设计, 各自输出)",
                        variable=self.mode, value="separate").pack(anchor="w", padx=6, pady=1)
        ttk.Radiobutton(mf, text="模式2 · 多文库共一张芯片(每文库一对专属引物, 单文件输出)",
                        variable=self.mode, value="shared").pack(anchor="w", padx=6, pady=1)

        right = ttk.Frame(pane)
        pane.add(right, weight=1)

        bf = ttk.Labelframe(right, text=" ② 基本参数 ", style="Step.TLabelframe")
        bf.pack(fill="x", padx=6, pady=6)
        self.e_flen = self._form_row(bf, 0, "F 长度 (nt)", "20", 6,
                                     "正向侧翼引物长度; 与 R 相同为单组输出")
        self.e_rlen = self._form_row(bf, 1, "R 长度 (nt)", "20", 6,
                                     "反向侧翼位点长度; 与 F 不同则分 F/R 两组")
        self.e_n = self._form_row(bf, 2, "每组生成条数", "20", 6,
                                  "每组的候选引物数量(模式2 固定每文库一对, 此项不适用)")
        self.e_seed = self._form_row(bf, 3, "随机种子 (可空)", "", 10,
                                     "填数字可精确复现同一结果; 留空每次随机")

        af = ttk.Labelframe(right, text=" 高级参数 (一般无需修改) ", style="Step.TLabelframe")
        af.pack(fill="x", padx=6, pady=6)
        self.e_gcmin = self._form_row(af, 0, "GC 含量下限", "0.40", 6)
        self.e_gcmax = self._form_row(af, 1, "GC 含量上限", "0.60", 6)
        self.e_dimer = self._form_row(af, 2, "二聚体检查长度", "8", 6,
                                      "引物间 3′ 端互补检查窗口; 0 = 关闭")
        self.e_clamp = self._form_row(af, 3, "3′ GC clamp (可空)", "", 6,
                                      "要求 3′ 端 5nt 内至少 N 个 G/C; 留空关闭")
        self.e_attempts = self._form_row(af, 4, "最大尝试次数", "2000000", 10,
                                         "每条引物的随机采样上限(防死循环)")
        self.e_thr = self._form_row(af, 5, "线程数", "4", 4,
                                    "并行进程; Windows 下自动退化为单线程")

        of = ttk.Labelframe(right, text=" 输出 (可选) ", style="Step.TLabelframe")
        of.pack(fill="x", padx=6, pady=6)
        self.e_out = self._form_path(of, 0, "输出文件/前缀", "",
                                     "单文件=引物文件名; 模式2=输出 tsv; 多文件模式1 留空自动命名")

        bar = ttk.Frame(self.tab_design)
        bar.pack(fill="x", padx=6, pady=6)
        self.btn_design = ttk.Button(bar, text="▶  开始设计", style="Primary.TButton",
                                     command=self.do_design)
        self.btn_design.pack(side="left", ipadx=14, ipady=4)
        ttk.Label(bar, style="Hint.TLabel",
                  text="完成后自动跳转『② 结果』页").pack(side="left", padx=10)

    def _form_row(self, parent, r, label, default, width, tip=None):
        ttk.Label(parent, text=label, width=14, anchor="e").grid(
            row=r, column=0, sticky="e", padx=(8, 4), pady=3)
        e = ttk.Entry(parent, width=width, font=FONT_MONO)
        e.insert(0, default)
        e.grid(row=r, column=1, sticky="w", padx=(0, 8), pady=3)
        if tip:
            Tooltip(e, tip)
        return e

    def _form_path(self, parent, r, label, default, tip=None, browse=None):
        ttk.Label(parent, text=label, width=14, anchor="e").grid(
            row=r, column=0, sticky="e", padx=(8, 4), pady=3)
        box = ttk.Frame(parent)
        box.grid(row=r, column=1, columnspan=2, sticky="we", padx=(0, 8), pady=3)
        e = ttk.Entry(box, font=FONT_MONO)
        e.insert(0, default)
        e.pack(side="left", fill="x", expand=True)
        if browse:
            ttk.Button(box, text="浏览…", width=8,
                       command=lambda: browse(e)).pack(side="left", padx=(4, 0))
        if tip:
            Tooltip(e, tip)
        return e


    # ---- 页2: 结果 ----
    def _build_results_tab(self):
        bar = ttk.Frame(self.tab_results)
        bar.pack(fill="x", padx=6, pady=(8, 2))
        self.res_hint = ttk.Label(bar, style="Hint.TLabel",
                                  text="尚无结果 — 请先在『① 设计引物』运行")
        self.res_hint.pack(side="left")
        ttk.Button(bar, text="复制选中", command=lambda: self._copy_tree_sel(self.res_tree)).pack(side="right", padx=2)
        ttk.Button(bar, text="复制全部", command=self.copy_results_all).pack(side="right", padx=2)
        ttk.Button(bar, text="打开输出文件夹", command=self.open_outdir).pack(side="right", padx=2)

        cols = ("name", "seq", "len", "gc", "tm")
        self.res_tree = ttk.Treeview(self.tab_results, columns=cols, show="headings",
                                     height=16, selectmode="extended")
        for cid, text, w, anchor in (("name", "名称", 70, "w"), ("seq", "序列", 300, "w"),
                                     ("len", "长度", 60, "e"), ("gc", "GC%", 60, "e"),
                                     ("tm", "Tm(℃)", 70, "e")):
            self.res_tree.heading(cid, text=text)
            self.res_tree.column(cid, width=w, anchor=anchor)
        self.res_tree.tag_configure("mono", font=FONT_MONO)
        ysb = ttk.Scrollbar(self.tab_results, orient="vertical", command=self.res_tree.yview)
        self.res_tree.configure(yscrollcommand=ysb.set)
        self.res_tree.pack(fill="both", expand=True, side="left", padx=(6, 0), pady=4)
        ysb.pack(fill="y", side="right", padx=(0, 6), pady=4)

        note = ttk.Label(self.tab_results, style="Hint.TLabel",
                         text="同长: 任取两条 P_i/P_j 作 F/R · 分长: F 组取正向, R 组取反向位点(R 的反向互补即反向引物) · 模式2: 每文库一行, 用各自 F/R 出合成单")
        note.pack(fill="x", padx=8, pady=(0, 6))

    # ---- 页3: 校验 ----
    def _build_check_tab(self):
        ff = ttk.Labelframe(self.tab_check, text=" 校验设置 ", style="Step.TLabelframe")
        ff.pack(fill="x", padx=6, pady=8)
        self.e_chk = self._form_path(ff, 0, "引物文件", "",
                                     "支持纯列表 / 设计输出 tsv / 模式2 shared tsv; 提取其中纯 ACGT 序列",
                                     self._browse_primer)
        ttk.Label(ff, style="Hint.TLabel",
                  text="校验对象 pool = 『① 设计引物』页文件列表中的全部文件(多文件按联合 pool, 即共芯片口径)"
                  ).grid(row=1, column=0, columnspan=3, sticky="w", padx=8, pady=(0, 6))

        bar = ttk.Frame(self.tab_check)
        bar.pack(fill="x", padx=6, pady=2)
        self.btn_check = ttk.Button(bar, text="▶  开始校验", style="Primary.TButton",
                                    command=self.do_check)
        self.btn_check.pack(side="left", ipadx=14, ipady=4)
        self.chk_hint = ttk.Label(bar, style="Hint.TLabel", text="")
        self.chk_hint.pack(side="left", padx=10)

        cols = ("name", "seq", "gc", "verdict", "reason")
        self.chk_tree = ttk.Treeview(self.tab_check, columns=cols, show="headings",
                                     height=12, selectmode="browse")
        for cid, text, w, anchor in (("name", "编号", 60, "w"), ("seq", "序列", 240, "w"),
                                     ("gc", "GC%", 55, "e"), ("verdict", "结果", 80, "center"),
                                     ("reason", "未通过原因", 460, "w")):
            self.chk_tree.heading(cid, text=text)
            self.chk_tree.column(cid, width=w, anchor=anchor)
        self.chk_tree.tag_configure("pass", foreground=C_PASS)
        self.chk_tree.tag_configure("fail", foreground=C_FAIL)
        self.chk_tree.tag_configure("mono", font=FONT_MONO)
        ysb = ttk.Scrollbar(self.tab_check, orient="vertical", command=self.chk_tree.yview)
        self.chk_tree.configure(yscrollcommand=ysb.set)
        self.chk_tree.pack(fill="both", expand=True, side="left", padx=(6, 0), pady=4)
        ysb.pack(fill="y", side="right", padx=(0, 6), pady=4)

    # ---- 页4: 合成序列 ----
    def _build_build_tab(self):
        ff = ttk.Labelframe(self.tab_build, text=" 生成设置 ", style="Step.TLabelframe")
        ff.pack(fill="x", padx=6, pady=8)
        self.e_bfasta = self._form_path(ff, 0, "payload FASTA", "",
                                        "待加侧翼的 oligo 序列", self._browse_fasta)
        self.e_pf = self._form_row(ff, 1, "引物 F (5′ 侧翼)", "", 26,
                                   "加在 5′ 端, 即 PCR 正向引物; 可从『② 结果』复制")
        self.e_pr = self._form_row(ff, 2, "引物 R (3′ 侧翼)", "", 26,
                                   "按原样加在 3′ 端; PCR 反向引物为其反向互补")
        self.rc_flag = tk.BooleanVar(value=False)
        ttk.Checkbutton(ff, text="给定的 R 已是反向引物本身 (自动取反向互补后拼接)",
                        variable=self.rc_flag).grid(row=3, column=1, sticky="w", padx=(0, 8), pady=2)
        self.e_bout = self._form_path(ff, 4, "输出前缀 (可空)", "",
                                      "默认在 fasta 同目录 oligo_for_synthesis")
        self.e_bmax = self._form_row(ff, 5, "长度告警 (nt)", "300", 6,
                                     "最终 oligo 超过该长度时提示(不截断)")

        bar = ttk.Frame(self.tab_build)
        bar.pack(fill="x", padx=6, pady=4)
        self.btn_build = ttk.Button(bar, text="▶  生成合成序列", style="Primary.TButton",
                                    command=self.do_build)
        self.btn_build.pack(side="left", ipadx=14, ipady=4)
        self.build_hint = ttk.Label(bar, style="Hint.TLabel", text="")
        self.build_hint.pack(side="left", padx=10)
        ttk.Button(bar, text="打开输出文件夹", command=self.open_outdir).pack(side="right", padx=2)
        ttk.Label(ff, style="Hint.TLabel",
                  text="生成前自动对当前 payload 校验两条引物(含反向互补形式); 不安全则拒绝写出"
                  ).grid(row=6, column=0, columnspan=3, sticky="w", padx=8, pady=(0, 6))

    # ---- 页5: 日志 ----
    def _build_log_tab(self):
        self.log_text = scrolledtext.ScrolledText(self.tab_log, font=(FONT_MONO[0], 9),
                                                  height=20, state="disabled", wrap="word")
        self.log_text.pack(fill="both", expand=True, padx=6, pady=6)
        bar = ttk.Frame(self.tab_log)
        bar.pack(fill="x", padx=6, pady=(0, 6))
        ttk.Button(bar, text="复制全部", command=self.copy_log).pack(side="left", padx=2)
        ttk.Button(bar, text="保存日志…", command=self.save_log).pack(side="left", padx=2)
        ttk.Button(bar, text="清空", command=lambda: self._log_set("")).pack(side="left", padx=2)


    # ---------- 文件浏览与通用 ----------
    def _browse_fasta(self, entry):
        p = filedialog.askopenfilename(filetypes=[("FASTA", "*.fa *.fasta *.fna"), ("全部", "*.*")])
        if p:
            entry.delete(0, "end")
            entry.insert(0, p)

    def _browse_primer(self, entry):
        p = filedialog.askopenfilename(title="选择引物文件",
                                       filetypes=[("引物/表格", "*.txt *.tsv *.csv"), ("全部", "*.*")])
        if p:
            entry.delete(0, "end")
            entry.insert(0, p)

    def add_files(self):
        ps = filedialog.askopenfilenames(
            filetypes=[("FASTA", "*.fa *.fasta *.fna"), ("全部", "*.*")])
        for p in ps:
            if p not in self.file_list.get(0, "end"):
                self.file_list.insert("end", p)

    def remove_sel(self):
        for i in reversed(self.file_list.curselection()):
            self.file_list.delete(i)

    def _log_set(self, s):
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.insert("end", s)
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def log(self, s):
        self.log_text.configure(state="normal")
        self.log_text.insert("end", s)
        self.log_text.see("end")
        self.log_text.configure(state="disabled")

    def copy_log(self):
        self.root.clipboard_clear()
        self.root.clipboard_append(self.log_text.get("1.0", "end"))

    def save_log(self):
        p = filedialog.asksaveasfilename(defaultextension=".txt",
                                         filetypes=[("文本", "*.txt")])
        if p:
            with open(p, "w", encoding="utf-8") as f:
                f.write(self.log_text.get("1.0", "end"))
            self.status.configure(text="日志已保存")

    def show_help(self):
        messagebox.showinfo("使用说明",
            "① 设计引物: 添加 FASTA(可多个) → 选模式 → 开始设计。k 阈值自动选择, 所有检查均含反向互补形式(实际反向引物)。" + chr(10) + chr(10) +
            "② 结果: 表格中挑选两条 — 同长任取 P_i/P_j; F/R 分长时 F 组取正向、R 组取反向位点; 模式2 每文库一行。" + chr(10) + chr(10) +
            "③ 校验: pool 更新后检查旧引物是否仍安全, 红/绿标示通过与否。" + chr(10) + chr(10) +
            "④ 生成合成序列: payload + F/R → fa/csv/pcr_primers.txt, 生成前自动复核安全性。" + chr(10) + chr(10) +
            "更多说明见 帮助 → 项目主页。")

    def show_about(self):
        messagebox.showinfo("关于",
            "%s v%s" % (APP_TITLE, APP_VERSION) + chr(10) + chr(10) +
            "基于 oligopool-primer-design 开源工具" + chr(10) + REPO_URL + chr(10) + chr(10) +
            "引物对所有 pool 序列(双链)做 k-mer 隔离, k 阈值按 pool 自动选择, 同一种子结果精确可复现。")

    def open_outdir(self):
        if not self._last_outputs:
            messagebox.showinfo(APP_TITLE, "还没有输出文件 — 请先运行设计或生成")
            return
        d = os.path.dirname(os.path.abspath(self._last_outputs[0]))
        try:
            os.startfile(d)
        except Exception as e:
            messagebox.showwarning(APP_TITLE, "无法打开: %r" % e)


    # ---------- 任务调度 ----------
    def _set_busy(self, busy, text):
        self.busy = busy
        for b in (self.btn_design, self.btn_check, self.btn_build):
            try:
                b.configure(state="disabled" if busy else "normal")
            except tk.TclError:
                pass
        if busy:
            self.progress.start(24)
            self.t_start = time.time()
        else:
            self.progress.stop()
            self.t_start = None
            self.elapsed.configure(text="")
        self.status.configure(text=text)

    def run_job(self, title, fn, on_done=None):
        if self.job is not None and self.job.is_alive():
            messagebox.showwarning(APP_TITLE, "有任务正在运行, 请稍候")
            return
        self.log("" + chr(10) + "===== %s =====" % title + chr(10))
        self._set_busy(True, "正在运行: %s" % title)
        self.job_done = on_done

        def worker():
            qw = QueueWriter(self.q)
            try:
                with contextlib.redirect_stdout(qw), contextlib.redirect_stderr(qw):
                    result = fn()
                self.q.put(None)
                self.q.put(("ok", result))
            except SystemExit as e:
                self.q.put(None)
                self.q.put(("exit", e.code))
            except Exception as e:
                self.q.put(None)
                self.q.put(("error", repr(e)))

        self.job = threading.Thread(target=worker, daemon=True)
        self.job.start()

    def _poll(self):
        try:
            while True:
                item = self.q.get_nowait()
                if item is None:
                    kind, payload = self.q.get_nowait()
                    self._finish_job(kind, payload)
                    continue
                if isinstance(item, tuple):
                    continue
                self.log(item)
        except queue.Empty:
            pass
        if self.t_start is not None:
            self.elapsed.configure(text="用时 %.0f s" % (time.time() - self.t_start))
        self.root.after(120, self._poll)

    def _finish_job(self, kind, payload):
        if kind == "ok":
            self.log(chr(10) + "[完成]" + chr(10))
            self._set_busy(False, "完成")
        elif kind == "exit":
            self.log(chr(10) + "[结束] 退出码 %s" % payload + chr(10))
            self._set_busy(False, "已结束 (退出码 %s)" % payload)
        else:
            self.log(chr(10) + "[出错] %s" % payload + chr(10))
            self._set_busy(False, "出错")
        cb, self.job_done = self.job_done, None
        if cb is not None:
            try:
                cb(kind)
            except Exception as e:
                self.log("[界面回调出错] %r" % e + chr(10))

    # ---------- 动作: 设计 ----------
    def _files(self):
        files = list(self.file_list.get(0, "end"))
        if not files:
            messagebox.showwarning(APP_TITLE, "请先在『① 设计引物』页添加 FASTA 文件")
            return None
        return files

    def do_design(self):
        files = self._files()
        if not files:
            return
        argv = ["--fasta"] + files + ["--mode", self.mode.get()]
        fl, rl = self.e_flen.get().strip(), self.e_rlen.get().strip()
        if fl and rl and fl != rl:
            argv += ["--primer-len", fl, rl]
        elif fl:
            argv += ["--primer-len", fl]
        for opt, ent in (("-n", self.e_n), ("--threads", self.e_thr),
                         ("--gc-min", self.e_gcmin), ("--gc-max", self.e_gcmax),
                         ("--dimer-k", self.e_dimer), ("--max-attempts", self.e_attempts)):
            v = ent.get().strip()
            if v:
                argv += [opt, v]
        if self.e_clamp.get().strip():
            argv += ["--gc-clamp", self.e_clamp.get().strip()]
        if self.e_seed.get().strip():
            argv += ["--seed", self.e_seed.get().strip()]
        if self.e_out.get().strip():
            argv += ["--out", self.e_out.get().strip()]
        shared = self.mode.get() == "shared"
        self.run_job("设计引物", lambda: opd.main(argv),
                     on_done=lambda kind: self._load_design_results(files, shared))

    def _load_design_results(self, files, shared):
        outs, rows = [], []
        if shared:
            out = self.e_out.get().strip() or os.path.join(
                os.path.dirname(os.path.abspath(files[0])), "safe_primers_shared.tsv")
            if os.path.exists(out):
                outs.append(out)
                for line in open(out, encoding="utf-8"):
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.rstrip(chr(10)).split(chr(9))
                    if parts and parts[0] in ("name", "library"):
                        continue
                    if len(parts) >= 4:
                        rows.append((parts[0], "F=" + parts[1] + "  R位=" + parts[2], "", "", ""))
        elif len(files) == 1:
            base = self.e_out.get().strip() or os.path.join(
                os.path.dirname(os.path.abspath(files[0])), "safe_primers.txt")
            tsv = base + ".tsv"
            if os.path.exists(tsv):
                outs.append(tsv)
                for line in open(tsv, encoding="utf-8"):
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.rstrip(chr(10)).split(chr(9))
                    if parts and parts[0] == "name":
                        continue
                    if len(parts) >= 5:
                        rows.append((parts[0], parts[1], parts[2], parts[3], parts[4]))
        else:
            for fp in files:
                stem = os.path.splitext(os.path.basename(fp))[0]
                tsv = os.path.join(os.path.dirname(os.path.abspath(fp)),
                                   stem + "_safe_primers.txt.tsv")
                if os.path.exists(tsv):
                    outs.append(tsv)
                    for line in open(tsv, encoding="utf-8"):
                        if line.startswith("#") or not line.strip():
                            continue
                        parts = line.rstrip(chr(10)).split(chr(9))
                        if parts and parts[0] == "name":
                            continue
                        if len(parts) >= 5:
                            rows.append((parts[0], parts[1], parts[2], parts[3], parts[4]))
        self._last_outputs = outs
        t = self.res_tree
        t.delete(*t.get_children())
        for r in rows:
            t.insert("", "end", values=r, tags=("mono",))
        if rows:
            self.nb.select(self.tab_results)
            self.res_hint.configure(
                text=("模式2 · 每文库一对 (共 %d 文库)" % len(rows)) if shared
                else ("共 %d 条候选 (来自 %d 个文件)" % (len(rows), len(files))))
        else:
            self.res_hint.configure(text="未找到输出文件 — 请查看日志页")

    def copy_results_all(self):
        t = self.res_tree
        lines = [chr(9).join(str(x) for x in t.item(iid, "values")) for iid in t.get_children("")]
        self.root.clipboard_clear()
        self.root.clipboard_append(chr(10).join(lines))
        self.status.configure(text="已复制 %d 行到剪贴板" % len(lines))

    def _copy_tree_sel(self, tree):
        sels = tree.selection()
        if not sels:
            return
        lines = [chr(9).join(str(x) for x in tree.item(iid, "values")) for iid in sels]
        self.root.clipboard_clear()
        self.root.clipboard_append(chr(10).join(lines))


    # ---------- 动作: 校验 ----------
    def do_check(self):
        files = self._files()
        if not files:
            return
        p = self.e_chk.get().strip()
        if not p or not os.path.exists(p):
            messagebox.showwarning(APP_TITLE, "请选择要校验的引物文件")
            return
        gcmin = float(self.e_gcmin.get() or 0.40)
        gcmax = float(self.e_gcmax.get() or 0.60)
        dk = int(self.e_dimer.get() or 8)

        def after(kind):
            self._fill_check()
            self.nb.select(self.tab_check)

        self.run_job("校验引物",
                     lambda: self._check_impl(files, p, gcmin, gcmax, dk),
                     on_done=after)

    def _check_impl(self, files, primer_path, gcmin, gcmax, dk):
        # 与命令行 --check 同一套判据, 直接产出结构化结果供表格(工作线程内运行)
        print("[读取] %d 个输入文件" % len(files))
        union = []
        for fp in files:
            seqs = [s.upper() for s in opd.read_fasta(fp) if s]
            print("  %s: %d 条" % (os.path.basename(fp), len(seqs)))
            union += seqs
        pool = union + [opd.rc_seq(s) for s in union]
        tk_, ik_, sets = opd.resolve_k_and_build(pool, primer_len=20)
        print("[判据] tail-k=%d / inner-k=%d (自动)" % (tk_, ik_))
        cfg = dict(gc_min=gcmin, gc_max=gcmax, max_homopolymer=5, max_gc_run=5,
                   gc_clamp=0, tail_k=tk_, inner_k=ik_, primer_len=20)
        primers, seen = [], set()
        for line in open(primer_path, encoding="utf-8"):
            for tok in line.split():
                s = tok.upper().strip(",;" + chr(34))
                if len(s) >= 10 and s not in seen and all(c in "ACGT" for c in s):
                    seen.add(s)
                    primers.append(s)
        self._check_rows = []
        n_fail = 0
        for i, p in enumerate(primers, 1):
            reasons = opd.fail_reasons_both(p, cfg, sets, detail=True)
            gc = 100.0 * (p.count("G") + p.count("C")) / len(p)
            verdict = "不通过" if reasons else "通过"
            if reasons:
                n_fail += 1
            self._check_rows.append(("P%02d" % i, p, "%.0f%%" % gc, verdict,
                                     "; ".join(reasons)))
        pair_fail = 0
        if dk > 0:
            for i in range(len(primers)):
                for j in range(i + 1, len(primers)):
                    if opd.rc_seq(primers[i]) == primers[j]:
                        continue
                    if opd.rc_seq(primers[i][-dk:]) in primers[j] or opd.rc_seq(primers[j][-dk:]) in primers[i]:
                        pair_fail += 1
                        self._check_rows.append(
                            ("配对", "P%02d + P%02d" % (i + 1, j + 1), "", "不通过",
                             "两条引物 3′ 端互补, 同管使用有二聚体风险"))
        print("[汇总] 通过 %d/%d; 二聚体风险配对 %d" % (len(primers) - n_fail, len(primers), pair_fail))
        if n_fail or pair_fail:
            raise SystemExit(1)

    def _fill_check(self):
        t = self.chk_tree
        t.delete(*t.get_children())
        for r in self._check_rows:
            tag = "fail" if r[3] == "不通过" else "pass"
            t.insert("", "end", values=r, tags=(tag, "mono"))
        n_pass = sum(1 for r in self._check_rows if r[3] == "通过")
        self.chk_hint.configure(text="通过 %d / %d 行" % (n_pass, len(self._check_rows)))

    # ---------- 动作: 生成 ----------
    def do_build(self):
        fa = self.e_bfasta.get().strip()
        pf = self.e_pf.get().strip().upper()
        pr = self.e_pr.get().strip().upper()
        if not (fa and pf and pr):
            messagebox.showwarning(APP_TITLE, "请填写 payload FASTA 与两条引物序列")
            return
        argv = ["--fasta", fa, "--primer-f", pf, "--primer-r", pr]
        if self.rc_flag.get():
            argv.append("--primer-r-is-rc")
        if self.e_bout.get().strip():
            argv += ["--out", self.e_bout.get().strip()]
        if self.e_bmax.get().strip():
            argv += ["--max-len", self.e_bmax.get().strip()]
        out = self.e_bout.get().strip() or os.path.join(
            os.path.dirname(os.path.abspath(fa)), "oligo_for_synthesis")
        self.run_job("生成合成序列", lambda: bo.main(argv),
                     on_done=lambda kind: self._after_build(out))

    def _after_build(self, out):
        made = [out + ext for ext in (".fa", ".csv", ".pcr_primers.txt") if os.path.exists(out + ext)]
        self._last_outputs = made
        if made:
            self.build_hint.configure(text="已生成: " + ", ".join(os.path.basename(x) for x in made))
        else:
            self.build_hint.configure(text="未生成输出 — 请查看日志页(校验未通过会拒绝写出)")


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
