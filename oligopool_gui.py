#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
oligopool_gui.py - Oligopool 引物设计工具图形界面

- 设计引物: 选 FASTA(可多个) + 模式 + 参数, 一键生成
- 校验引物: 对当前 pool 检查已有引物文件是否仍安全
- 生成合成序列: payload + F/R 引物 -> fa / csv / pcr_primers.txt
依赖同目录的 oligopool_primer_design.py 与 build_oligos.py
"""
import contextlib
import queue
import threading
import tkinter as tk
from tkinter import ttk, filedialog, scrolledtext, messagebox

import oligopool_primer_design as opd
import build_oligos as bo

APP_TITLE = "Oligopool 引物设计工具"


class QueueWriter(object):
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
        root.title(APP_TITLE)
        root.geometry("880x760")
        self.q = queue.Queue()
        self.job = None
        self._build()
        self.root.after(100, self._poll)

    # ---------- 界面 ----------
    def _build(self):
        top = ttk.Frame(self.root)
        top.pack(fill="both", expand=True)

        ff = ttk.LabelFrame(top, text="输入 FASTA 文件(可多选; 多文件时按所选模式处理)")
        ff.pack(fill="x", padx=8, pady=4)
        self.file_list = tk.Listbox(ff, height=4, selectmode="extended")
        self.file_list.pack(side="left", fill="both", expand=True,
                            padx=(6, 2), pady=4)
        fb = ttk.Frame(ff)
        fb.pack(side="left", fill="y", padx=4)
        ttk.Button(fb, text="添加文件...", command=self.add_files).pack(fill="x", pady=2)
        ttk.Button(fb, text="移除所选", command=self.remove_sel).pack(fill="x", pady=2)

        mf = ttk.LabelFrame(top, text="多文件模式")
        mf.pack(fill="x", padx=8, pady=4)
        self.mode = tk.StringVar(value="separate")
        ttk.Radiobutton(mf, text="模式1: 每个文件单独一张芯片, 依次独立设计, 各自输出",
                        variable=self.mode, value="separate").pack(anchor="w", padx=6)
        ttk.Radiobutton(mf, text="模式2: 多文库共一张芯片, 每文库一对专属引物, 输出单个文件",
                        variable=self.mode, value="shared").pack(anchor="w", padx=6)

        pf = ttk.LabelFrame(top, text="设计参数")
        pf.pack(fill="x", padx=8, pady=4)
        r1 = ttk.Frame(pf)
        r1.pack(fill="x", pady=3)
        self.e_flen = self._entry(r1, "F 长度(nt)", "20", 5)
        self.e_rlen = self._entry(r1, "R 长度(nt)", "20", 5)
        self.e_n = self._entry(r1, "每组条数", "20", 5)
        self.e_seed = self._entry(r1, "随机种子(可空)", "", 8)
        self.e_thr = self._entry(r1, "线程", "4", 4)
        r2 = ttk.Frame(pf)
        r2.pack(fill="x", pady=3)
        self.e_gcmin = self._entry(r2, "GC 下限", "0.40", 6)
        self.e_gcmax = self._entry(r2, "GC 上限", "0.60", 6)

        ab = ttk.Frame(top)
        ab.pack(fill="x", padx=8, pady=4)
        ttk.Button(ab, text="★ 设计引物", command=self.do_design).pack(
            side="left", padx=4, ipadx=10, ipady=3)
        ttk.Button(ab, text="校验已有引物...", command=self.do_check).pack(
            side="left", padx=4, ipady=3)

        bf = ttk.LabelFrame(top, text="生成送合成的完整序列 (build_oligos)")
        bf.pack(fill="x", padx=8, pady=4)
        r3 = ttk.Frame(bf)
        r3.pack(fill="x", pady=2)
        self.e_bfasta = self._entry(r3, "payload FASTA", "", 34)
        ttk.Button(r3, text="浏览...", command=lambda: self._browse(self.e_bfasta)).pack(
            side="left", padx=2)
        r4 = ttk.Frame(bf)
        r4.pack(fill="x", pady=2)
        self.e_pf = self._entry(r4, "引物F(5′侧翼)", "", 24)
        self.e_pr = self._entry(r4, "引物R(3′侧翼)", "", 24)
        self.rc_flag = tk.BooleanVar(value=False)
        ttk.Checkbutton(r4, text="R 已给成反向引物(自动取反补)",
                        variable=self.rc_flag).pack(side="left", padx=4)
        r5 = ttk.Frame(bf)
        r5.pack(fill="x", pady=2)
        self.e_bout = self._entry(r5, "输出前缀(可空)", "", 30)
        ttk.Button(bf, text="生成合成序列", command=self.do_build).pack(
            anchor="w", padx=6, pady=3, ipady=2)

        lf = ttk.LabelFrame(top, text="运行日志")
        lf.pack(fill="both", expand=True, padx=8, pady=4)
        self.log_text = scrolledtext.ScrolledText(lf, height=14, font=("Consolas", 9))
        self.log_text.pack(fill="both", expand=True, padx=4, pady=4)
        self.log("使用说明: 先添加 FASTA -> 选模式 -> [设计引物];"
                 " 之后把选定的 F/R 填到下方生成合成序列。\n"
                 "k 阈值(tail-k/inner-k)自动选择; 全部检查均含反向互补形式。\n")

    def _entry(self, parent, label, default, width):
        f = ttk.Frame(parent)
        f.pack(side="left", padx=6)
        ttk.Label(f, text=label).pack(side="top", anchor="w")
        e = ttk.Entry(f, width=width)
        e.insert(0, default)
        e.pack(side="top")
        return e

    def _browse(self, entry):
        p = filedialog.askopenfilename(filetypes=[("FASTA", "*.fa *.fasta *.fna"), ("全部", "*.*")])
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

    def log(self, s):
        self.log_text.insert("end", s)
        self.log_text.see("end")

    # ---------- 任务调度 ----------
    def run_job(self, title, fn):
        if self.job is not None and self.job.is_alive():
            messagebox.showwarning(APP_TITLE, "有任务正在运行, 请稍候")
            return
        self.log("\n===== %s =====\n" % title)

        def worker():
            qw = QueueWriter(self.q)
            try:
                with contextlib.redirect_stdout(qw), contextlib.redirect_stderr(qw):
                    fn()
                self.q.put("\n[完成] %s\n" % title)
            except SystemExit as e:
                self.q.put("\n[结束] %s (退出码 %s)\n" % (title, e.code))
            except Exception as e:
                self.q.put("\n[出错] %r\n" % (e,))

        self.job = threading.Thread(target=worker, daemon=True)
        self.job.start()

    def _poll(self):
        try:
            while True:
                self.log(self.q.get_nowait())
        except queue.Empty:
            pass
        self.root.after(120, self._poll)

    # ---------- 动作 ----------
    def _files(self):
        files = list(self.file_list.get(0, "end"))
        if not files:
            messagebox.showwarning(APP_TITLE, "请先添加 FASTA 文件")
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
        if self.e_n.get().strip():
            argv += ["-n", self.e_n.get().strip()]
        if self.e_seed.get().strip():
            argv += ["--seed", self.e_seed.get().strip()]
        if self.e_thr.get().strip():
            argv += ["--threads", self.e_thr.get().strip()]
        if self.e_gcmin.get().strip() and self.e_gcmax.get().strip():
            argv += ["--gc-min", self.e_gcmin.get().strip(),
                     "--gc-max", self.e_gcmax.get().strip()]
        self.run_job("设计引物 %s" % ("; ".join(files)), lambda: opd.main(argv))

    def do_check(self):
        files = self._files()
        if not files:
            return
        p = filedialog.askopenfilename(title="选择要校验的引物文件",
                                       filetypes=[("文本", "*.txt *.tsv *.csv"), ("全部", "*.*")])
        if not p:
            return
        argv = ["--fasta"] + files + ["--check", p]
        self.run_job("校验引物 %s" % p, lambda: opd.main(argv))

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
        self.run_job("生成合成序列 %s" % fa, lambda: bo.main(argv))


def main():
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
