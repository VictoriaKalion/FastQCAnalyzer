# -*- coding: utf-8 -*-
import subprocess
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from PIL import Image, ImageTk, ImageDraw
import matplotlib.pyplot as plt
from collections import defaultdict
import os

print("УСТАНАВЛИВАЕМ MATPLOTLIB...")
subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib"])

print("ЗАПУСКАЕМ FASTQ АНАЛИЗАТОР...")
import matplotlib.pyplot as plt
from collections import defaultdict


class FastqReader:
    """
    Класс для чтения и анализа FASTQ файлов с оптимизацией памяти
    Использует генераторы для работы с большими файлами
    """

    def __init__(self, filename):
        self.filename = filename
        self._sequence_count = None
        self._total_length = None

    def _read_fastq_chunks(self):
        """ГЕНЕРАТОР: читает FASTQ файл по одному риду за раз"""
        with open(self.filename, 'r', encoding='utf-8') as file:
            while True:
                # Читаем 4 строки одного рида
                lines = [file.readline().strip() for _ in range(4)]
                if not lines[0]:  # Конец файла
                    break
                yield lines

    def calculate_statistics(self):
        """Рассчитывает статистику используя генератор (память O(1))"""
        count = 0
        total_length = 0

        for chunk in self._read_fastq_chunks():
            count += 1
            sequence_line = chunk[1]  # Вторая строка - последовательность
            total_length += len(sequence_line)

        self._sequence_count = count
        self._total_length = total_length
        return count, total_length

    def get_sequence_count(self):
        """Возвращает количество последовательностей в файле"""
        if self._sequence_count is None:
            self.calculate_statistics()
        return self._sequence_count

    def get_average_length(self):
        """Возвращает среднюю длину последовательностей"""
        if self._sequence_count is None:
            self.calculate_statistics()
        if self._sequence_count == 0:
            return 0
        return self._total_length / self._sequence_count

    def plot_per_base_quality(self, output="quality.png"):
        """Строит график качества по позициям (Per Base Sequence Quality)"""
        quality_sums = defaultdict(int)
        quality_counts = defaultdict(int)
        max_position = 0

        for chunk in self._read_fastq_chunks():
            quality_line = chunk[3]  # Четвертая строка - качество
            for i, char in enumerate(quality_line):
                score = ord(char) - 33  # Конвертируем в числовое качество
                quality_sums[i] += score
                quality_counts[i] += 1
                max_position = max(max_position, i)

        # Рассчитываем среднее качество для каждой позиции
        positions = range(1, max_position + 2)
        avg_qualities = [quality_sums[i] / quality_counts[i] for i in range(max_position + 1)]

        # Создаем график в монохромной гамме
        plt.figure(figsize=(10, 6), facecolor='#FFFFFF')
        ax = plt.gca()
        ax.set_facecolor('#FAFAFA')

        plt.plot(positions, avg_qualities, linewidth=3, color='#424242', marker='o',
                 markersize=5, markerfacecolor='#212121', markeredgecolor='#212121',
                 markeredgewidth=1)

        plt.title('Качество последовательностей по позициям', fontfamily='Georgia',
                  color='#212121', fontsize=14, pad=20, style='italic')
        plt.xlabel('Позиция в риде (bp)', fontfamily='Georgia', color='#424242', fontsize=11)
        plt.ylabel('Phred Quality Score', fontfamily='Georgia', color='#424242', fontsize=11)

        # Настраиваем цвета осей и сетки
        ax.tick_params(colors='#616161', which='both')
        ax.spines['bottom'].set_color('#BDBDBD')
        ax.spines['top'].set_color('#BDBDBD')
        ax.spines['right'].set_color('#BDBDBD')
        ax.spines['left'].set_color('#BDBDBD')
        ax.grid(True, alpha=0.4, color='#E0E0E0', linestyle='--')

        plt.tight_layout()
        plt.savefig(output, dpi=150, facecolor='#FFFFFF', edgecolor='none', bbox_inches='tight')
        plt.close()
        return output

    def plot_per_base_content(self, output="content.png"):
        """Строит график содержания нуклеотидов по позициям (Per Base Sequence Content)"""
        base_counts = {'A': defaultdict(int), 'C': defaultdict(int),
                       'G': defaultdict(int), 'T': defaultdict(int)}
        total_counts = defaultdict(int)
        max_position = 0

        for chunk in self._read_fastq_chunks():
            sequence_line = chunk[1]  # Вторая строка - последовательность
            for i, base in enumerate(sequence_line):
                base_upper = base.upper()
                if base_upper in base_counts:
                    base_counts[base_upper][i] += 1
                    total_counts[i] += 1
                    max_position = max(max_position, i)

        positions = range(1, max_position + 2)
        plt.figure(figsize=(10, 6), facecolor='#FFFFFF')
        ax = plt.gca()
        ax.set_facecolor('#FAFAFA')

        # Монохромная палитра для нуклеотидов
        colors = ['#212121', '#424242', '#616161', '#757575', '#9E9E9E', '#BDBDBD']

        # Строим линии для каждого нуклеотида
        for (base, counts), color in zip(base_counts.items(), colors):
            percentages = [counts[i] / total_counts[i] * 100 if total_counts[i] > 0 else 0
                           for i in range(max_position + 1)]
            plt.plot(positions, percentages, label=base, linewidth=2.5, color=color)

        plt.title('Содержание нуклеотидов по позициям', fontfamily='Georgia',
                  color='#212121', fontsize=14, pad=20, style='italic')
        plt.xlabel('Позиция в риде (bp)', fontfamily='Georgia', color='#424242', fontsize=11)
        plt.ylabel('Процент (%)', fontfamily='Georgia', color='#424242', fontsize=11)

        # Легенда с монохромными цветами
        legend = plt.legend(prop={'family': 'Georgia', 'size': 10})
        for text in legend.get_texts():
            text.set_color('#424242')

        # Настраиваем цвета осей и сетки
        ax.tick_params(colors='#616161', which='both')
        ax.spines['bottom'].set_color('#BDBDBD')
        ax.spines['top'].set_color('#BDBDBD')
        ax.spines['right'].set_color('#BDBDBD')
        ax.spines['left'].set_color('#BDBDBD')
        ax.grid(True, alpha=0.4, color='#E0E0E0', linestyle='--')

        plt.tight_layout()
        plt.savefig(output, dpi=150, facecolor='#FFFFFF', edgecolor='none', bbox_inches='tight')
        plt.close()
        return output

    def plot_sequence_length_distribution(self, output="length.png"):
        """Строит гистограмму распределения длин последовательностей"""
        lengths = []

        for chunk in self._read_fastq_chunks():
            sequence_line = chunk[1]  # Вторая строка - последовательность
            lengths.append(len(sequence_line))

        plt.figure(figsize=(10, 6), facecolor='#FFFFFF')
        ax = plt.gca()
        ax.set_facecolor('#FAFAFA')

        # Монохромные цвета для гистограммы
        plt.hist(lengths, bins=20, edgecolor='#424242', alpha=0.8, color='#757575')

        plt.title('Распределение длин последовательностей', fontfamily='Georgia',
                  color='#212121', fontsize=14, pad=20, style='italic')
        plt.xlabel('Длина последовательности (bp)', fontfamily='Georgia', color='#424242', fontsize=11)
        plt.ylabel('Частота', fontfamily='Georgia', color='#424242', fontsize=11)

        # Настраиваем цвета осей и сетки
        ax.tick_params(colors='#616161', which='both')
        ax.spines['bottom'].set_color('#BDBDBD')
        ax.spines['top'].set_color('#BDBDBD')
        ax.spines['right'].set_color('#BDBDBD')
        ax.spines['left'].set_color('#BDBDBD')
        ax.grid(True, alpha=0.4, color='#E0E0E0', linestyle='--')

        plt.tight_layout()
        plt.savefig(output, dpi=150, facecolor='#FFFFFF', edgecolor='none', bbox_inches='tight')
        plt.close()
        return output


class RoundedButton(tk.Canvas):
    """Кастомная скругленная кнопка на основе Canvas"""

    def __init__(self, parent, text, command, width=150, height=45,
                 corner_radius=20, bg_color='#E0E0E0', hover_color='#BDBDBD',
                 text_color='#212121', font=('Georgia', 10, 'bold')):
        super().__init__(parent, width=width, height=height,
                         highlightthickness=0, bg=parent.cget('bg'))

        self.command = command
        self.bg_color = bg_color
        self.hover_color = hover_color
        self.text_color = text_color
        self.font = font
        self.corner_radius = corner_radius
        self.width = width
        self.height = height
        self.text = text

        self.bind("<Button-1>", self._on_click)
        self.bind("<Enter>", self._on_enter)
        self.bind("<Leave>", self._on_leave)

        self.draw_button(bg_color)

    def draw_button(self, color):
        """Рисует скругленную кнопку"""
        self.delete("all")

        # Рисуем скругленный прямоугольник с увеличенным радиусом
        self.create_rounded_rect(2, 2, self.width - 2, self.height - 2,
                                 self.corner_radius, fill=color, outline='#BDBDBD', width=1)

        # Добавляем текст
        self.create_text(self.width // 2, self.height // 2,
                         text=self.text, fill=self.text_color, font=self.font)

    def create_rounded_rect(self, x1, y1, x2, y2, radius, **kwargs):
        """Создает скругленный прямоугольник"""
        points = [x1 + radius, y1,
                  x2 - radius, y1,
                  x2, y1,
                  x2, y1 + radius,
                  x2, y2 - radius,
                  x2, y2,
                  x2 - radius, y2,
                  x1 + radius, y2,
                  x1, y2,
                  x1, y2 - radius,
                  x1, y1 + radius,
                  x1, y1]

        return self.create_polygon(points, smooth=True, **kwargs)

    def _on_click(self, event):
        self.command()

    def _on_enter(self, event):
        self.draw_button(self.hover_color)

    def _on_leave(self, event):
        self.draw_button(self.bg_color)


class FastQCAnalyzerGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("FastQC Analyzer")
        self.root.geometry("900x650")
        self.root.configure(bg='#FFFFFF')

        # Установка шрифта Georgia с курсивом для заголовков
        self.style = ttk.Style()
        self.style.configure('TLabel', font=('Georgia', 11), background='#FFFFFF', foreground='#212121')
        self.style.configure('TButton', font=('Georgia', 10), borderwidth=0, focuscolor='none')
        self.style.configure('Title.TLabel', font=('Georgia', 16, 'bold', 'italic'))
        self.style.configure('Italic.TLabel', font=('Georgia', 11, 'italic'))

        self.filename = None
        self.analyzer = None

        self.setup_ui()

    def setup_ui(self):
        # Главный фрейм
        main_frame = tk.Frame(self.root, bg='#FFFFFF')
        main_frame.pack(fill='both', expand=True, padx=25, pady=25)

        # Заголовок
        title_label = ttk.Label(main_frame, text="FastQC Analyzer", style='Title.TLabel')
        title_label.pack(pady=(0, 25))

        # Фрейм для выбора файла
        file_frame = tk.Frame(main_frame, bg='#FAFAFA', bd=1, relief='solid',
                              highlightbackground='#E0E0E0', highlightthickness=1)
        file_frame.pack(fill='x', pady=12, padx=10, ipady=12)

        ttk.Label(file_frame, text="Выберите FASTQ файл:", style='Italic.TLabel').pack(side='left', padx=15)

        # Скругленная кнопка выбора файла с увеличенным радиусом
        self.select_btn = RoundedButton(
            file_frame, "Выбрать файл", self.select_file,
            width=150, height=45, corner_radius=20,
            bg_color='#F5F5F5', hover_color='#E0E0E0', text_color='#212121'
        )
        self.select_btn.pack(side='left', padx=15)

        # Метка с именем файла
        self.file_label = ttk.Label(file_frame, text="Файл не выбран", foreground='#616161', style='Italic.TLabel')
        self.file_label.pack(side='left', padx=15, fill='x', expand=True)

        # Фрейм статистики
        stats_frame = tk.Frame(main_frame, bg='#FAFAFA', bd=1, relief='solid',
                               highlightbackground='#E0E0E0', highlightthickness=1)
        stats_frame.pack(fill='x', pady=12, padx=10, ipady=12)

        ttk.Label(stats_frame, text="Статистика файла:", style='Title.TLabel').pack(anchor='w', padx=15, pady=8)

        self.stats_text = tk.Text(stats_frame, height=6, width=80, font=('Georgia', 10),
                                  bg='#FFFFFF', fg='#212121', relief='flat', wrap='word',
                                  borderwidth=1, highlightthickness=1, highlightcolor='#E0E0E0',
                                  selectbackground='#E0E0E0')
        self.stats_text.pack(padx=15, pady=10, fill='x')
        self.stats_text.insert('1.0', "Статистика будет отображена здесь после выбора файла...")
        self.stats_text.config(state='disabled')

        # Фрейм для кнопок графиков
        plots_frame = tk.Frame(main_frame, bg='#FAFAFA', bd=1, relief='solid',
                               highlightbackground='#E0E0E0', highlightthickness=1)
        plots_frame.pack(fill='x', pady=12, padx=10, ipady=12)

        ttk.Label(plots_frame, text="Генерация графиков:", style='Title.TLabel').pack(anchor='w', padx=15, pady=8)

        buttons_frame = tk.Frame(plots_frame, bg='#FAFAFA')
        buttons_frame.pack(padx=15, pady=12, fill='x')

        # Монохромная палитра для кнопок
        button_colors = [
            ('#F5F5F5', '#E0E0E0'),  # Светло-серый -> Серый
            ('#EEEEEE', '#E0E0E0'),  # Очень светло-серый -> Серый
            ('#FAFAFA', '#E0E0E0'),  # Бежево-белый -> Серый
            ('#E0E0E0', '#BDBDBD')  # Серый -> Темно-серый
        ]

        # Скругленные кнопки для генерации графиков с увеличенным радиусом
        self.quality_btn = RoundedButton(
            buttons_frame, "Качество", self.plot_quality,
            width=150, height=45, corner_radius=20,
            bg_color=button_colors[0][0], hover_color=button_colors[0][1], text_color='#212121'
        )
        self.quality_btn.pack(side='left', padx=8, pady=6)

        self.content_btn = RoundedButton(
            buttons_frame, "Нуклеотиды", self.plot_content,
            width=150, height=45, corner_radius=20,
            bg_color=button_colors[1][0], hover_color=button_colors[1][1], text_color='#212121'
        )
        self.content_btn.pack(side='left', padx=8, pady=6)

        self.length_btn = RoundedButton(
            buttons_frame, "Длины", self.plot_length,
            width=150, height=45, corner_radius=20,
            bg_color=button_colors[2][0], hover_color=button_colors[2][1], text_color='#212121'
        )
        self.length_btn.pack(side='left', padx=8, pady=6)

        self.all_plots_btn = RoundedButton(
            buttons_frame, "Все графики", self.plot_all,
            width=150, height=45, corner_radius=20,
            bg_color=button_colors[3][0], hover_color=button_colors[3][1], text_color='#212121'
        )
        self.all_plots_btn.pack(side='left', padx=8, pady=6)

        # Область для вывода изображений
        self.image_frame = tk.Frame(main_frame, bg='#FFFFFF')
        self.image_frame.pack(fill='both', expand=True, pady=15)

        self.image_label = ttk.Label(self.image_frame, text="Графики будут отображены здесь",
                                     background='#FFFFFF', foreground='#757575', style='Italic.TLabel')
        self.image_label.pack(expand=True)

        # Блокировка кнопок до выбора файла
        self.set_buttons_state('disabled')

    def set_buttons_state(self, state):
        """Устанавливает состояние кнопок"""
        for btn in [self.quality_btn, self.content_btn, self.length_btn, self.all_plots_btn]:
            if state == 'disabled':
                # Визуально делаем кнопки светлыми когда отключены
                btn.bg_color = '#FAFAFA'
                btn.hover_color = '#FAFAFA'
                btn.draw_button('#FAFAFA')
            else:
                # Возвращаем оригинальные цвета
                if btn == self.quality_btn:
                    btn.bg_color = '#F5F5F5'
                    btn.hover_color = '#E0E0E0'
                elif btn == self.content_btn:
                    btn.bg_color = '#EEEEEE'
                    btn.hover_color = '#E0E0E0'
                elif btn == self.length_btn:
                    btn.bg_color = '#FAFAFA'
                    btn.hover_color = '#E0E0E0'
                elif btn == self.all_plots_btn:
                    btn.bg_color = '#E0E0E0'
                    btn.hover_color = '#BDBDBD'
                btn.draw_button(btn.bg_color)

    def select_file(self):
        """Выбор файла через диалоговое окно"""
        filename = filedialog.askopenfilename(
            title="Выберите FASTQ файл",
            filetypes=[("FASTQ files", "*.fastq *.fq"), ("All files", "*.*")]
        )

        if filename:
            self.filename = filename
            self.file_label.config(text=os.path.basename(filename))
            self.analyze_file()

    def analyze_file(self):
        """Анализ выбранного файла"""
        try:
            self.analyzer = FastqReader(self.filename)
            count = self.analyzer.get_sequence_count()
            avg_len = self.analyzer.get_average_length()
            total_bp = count * avg_len

            stats_text = f"""Количество последовательностей: {count:,}
Средняя длина: {avg_len:.2f} bp
Общий объем данных: {total_bp:,.0f} bp
Размер файла: {os.path.getsize(self.filename) / 1024 / 1024:.2f} MB"""

            self.stats_text.config(state='normal')
            self.stats_text.delete('1.0', 'end')
            self.stats_text.insert('1.0', stats_text)
            self.stats_text.config(state='disabled')

            self.set_buttons_state('normal')

        except Exception as e:
            messagebox.showerror("Ошибка", f"Ошибка при анализе файла: {str(e)}")

    def plot_quality(self):
        """Генерация графика качества"""
        if self.analyzer:
            try:
                output_file = self.analyzer.plot_per_base_quality()
                self.display_image(output_file)
            except Exception as e:
                messagebox.showerror("Ошибка", f"Ошибка при создании графика: {str(e)}")

    def plot_content(self):
        """Генерация графика содержания нуклеотидов"""
        if self.analyzer:
            try:
                output_file = self.analyzer.plot_per_base_content()
                self.display_image(output_file)
            except Exception as e:
                messagebox.showerror("Ошибка", f"Ошибка при создании графика: {str(e)}")

    def plot_length(self):
        """Генерация графика распределения длин"""
        if self.analyzer:
            try:
                output_file = self.analyzer.plot_sequence_length_distribution()
                self.display_image(output_file)
            except Exception as e:
                messagebox.showerror("Ошибка", f"Ошибка при создании графика: {str(e)}")

    def plot_all(self):
        """Генерация всех графиков"""
        if self.analyzer:
            try:
                quality_file = self.analyzer.plot_per_base_quality()
                content_file = self.analyzer.plot_per_base_content()
                length_file = self.analyzer.plot_sequence_length_distribution()

                # Показываем последний созданный график
                self.display_image(length_file)
                messagebox.showinfo("Успех", "Все графики созданы успешно!")

            except Exception as e:
                messagebox.showerror("Ошибка", f"Ошибка при создании графиков: {str(e)}")

    def display_image(self, image_path):
        """Отображение изображения в интерфейсе"""
        try:
            image = Image.open(image_path)
            # Масштабируем изображение под размер окна
            width, height = image.size
            max_width = 800
            max_height = 400

            if width > max_width or height > max_height:
                ratio = min(max_width / width, max_height / height)
                new_size = (int(width * ratio), int(height * ratio))
                image = image.resize(new_size, Image.Resampling.LANCZOS)

            photo = ImageTk.PhotoImage(image)

            # Очищаем предыдущее изображение
            for widget in self.image_frame.winfo_children():
                widget.destroy()

            # Создаем новую метку с изображением
            img_label = tk.Label(self.image_frame, image=photo, bg='#FFFFFF')
            img_label.image = photo  # Сохраняем ссылку
            img_label.pack(expand=True)

        except Exception as e:
            messagebox.showerror("Ошибка", f"Не удалось загрузить изображение: {str(e)}")


def create_test_fastq():
    """Создает тестовый FASTQ файл для демонстрации"""
    print("СОЗДАЕМ ТЕСТОВЫЙ ФАЙЛ...")
    test_data = ""
    for i in range(100):
        test_data += f"@read{i}\n"
        test_data += "ATCGATCGATCGATCGATCG\n"  # 20 bp
        test_data += "+\n"
        test_data += "IIIIIIIIIIIIIIIIIIII\n"  # Качество 40
    return test_data


if __name__ == "__main__":
    # Создаем тестовый файл если нужно
    test_data = create_test_fastq()
    with open("test.fastq", "w", encoding='utf-8') as f:
        f.write(test_data)

    # Запускаем GUI
    root = tk.Tk()
    app = FastQCAnalyzerGUI(root)
    root.mainloop()