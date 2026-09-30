# Script for FCS-GX decontaminator result visualisation

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

action_report = "/home/kiselevm/Рабочий стол/abyss-exclude.txt"
taxonomy_report = "/home/kiselevm/Рабочий стол/abyss-taxonomy-report.txt"
assembler='ABySS'

# Настройка шрифтов для корректного отображения русского языка в графиках
plt.rcParams['font.family'] = 'sans-serif'
# Для Windows обычно подходит 'Arial', для Mac - 'Arial', для Linux - 'DejaVu Sans'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']
plt.rcParams['axes.unicode_minus'] = False

# To rename those on the plot:
taxa_dict = {
    'anml':'Animals',
    'plnt':'Plants',
    'prok':'Prokaryotes',
    'prst':'Protista',
    'virs':'Viruses'
}

def plot_fcs_coverage(action_report, taxonomy_report, title):
    # 1. Загружаем данные

    df = pd.read_csv(action_report, sep='\t')
    df_tax = pd.read_csv(taxonomy_report, sep='\t', dtype=str, skiprows=1)


    # 2. Извлекаем главную группу из колонки 'div' (до двоеточия)
    df['div_group'] = df['div'].str.split(':').str[0].replace(taxa_dict)

    # 3. Убеждаемся, что agg_cont_cov числовой (на всякий случай отбрасываем возможные 'n/a')
    df['agg_cont_cov'] = pd.to_numeric(df['agg_cont_cov'], errors='coerce')
    df = df.dropna(subset=['agg_cont_cov'])

    # 4. Разбиваем значения на диапазоны (бинируем)
    # Правая граница включается, поэтому интервалы: (0, 50], (50, 90], (90, 100]
    bins = [0, 50, 90, 100]
    labels = ['0-50%', '50-90%', '90-100%']
    df['cov_range'] = pd.cut(df['agg_cont_cov'], bins=bins, labels=labels, right=True, include_lowest=True)

    # Добавляем информацию по alignment score
    # (берём только первое выравнивание из 4-х, самое первое)
    # (пока что это при построении графика не используется)
    score_map = df_tax.set_index('seq-id')['score-1'].dropna().to_dict()
    df['alignment_score'] = df['seq-id'].astype(str).map(score_map)

    # 5. Группируем данные: считаем количество записей для каждой группы в каждом диапазоне
    # unstack() превращает таблицу в удобный формат для столбчатой диаграммы
    counts = df.groupby(['div_group', 'cov_range']).size().unstack(fill_value=0)

    # Убедимся, что все три колонки диапазонов присутствуют, даже если в них нули
    for label in labels:
        if label not in counts.columns:
            counts[label] = 0

    # Сортируем колонки в правильном порядке
    counts = counts[labels]

    # 6. Построение графика
    fig, ax = plt.subplots(figsize=(10, 6))

    # Строим столбчатую диаграмму с группировкой
    counts.plot(kind='bar', ax=ax, width=0.8, color=['#d9534f', '#f0ad4e', '#5cb85c'], edgecolor='black')

    # Настройки внешнего вида
    ax.set_title(title, fontsize=14, pad=15)
    ax.set_xlabel('Таксономическая группа', fontsize=12)
    ax.set_ylabel('Количество последовательностей', fontsize=12)

    # Поворачиваем подписи оси X для читаемости
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, ha='center')

    # Настройка легенды
    ax.legend(title='Процент покрытия\nвыравниванием последовательности', bbox_to_anchor=(1.02, 1), loc='upper left')

    # Добавляем сетку для легкого чтения значений по оси Y
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)

    # Убираем верхнюю и правую рамки графика для стиля
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Подгоняем отступы, чтобы легенда не обрезалась
    plt.tight_layout()

    # Сохраняем график в файл и показываем на экране
    plt.savefig(f'{title}-coverage.png', dpi=300, bbox_inches='tight')

    ##########################################################################

    # LENGTH DISTRIBUTION
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_title(title)
    # Строим boxplot.
    # showfliers=True показывает выбросы (точки за "усами")
    sns.boxplot(
        data=df,
        x='div_group',
        y='seq_len',
        palette='pastel',
        showfliers=True,
        ax=ax
    )
    ax.set_xlabel('Таксономическая группа', fontsize=12)
    ax.set_ylabel('Длина последовательности (п.н.)', fontsize=12)

    # Сетка для легкого чтения
    ax.yaxis.grid(True, linestyle='--', alpha=0.7)
    ax.set_axisbelow(True)

    # Убираем рамки
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_ylim(0, 25000)

    plt.tight_layout()
    plt.savefig(f'{title}-length.png', dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    action_report = action_report
    taxonomy_report = taxonomy_report
    assembler=assembler
    plot_fcs_coverage(action_report, taxonomy_report, assembler)