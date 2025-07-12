

# Импорты: стандартные, сторонние, локальные

import os
import streamlit as st
import asyncio
import httpx
import plotly.express as px
from app.parsers import load_genotype_from_txt, load_genotype_from_vcf
from app.logic import annotate_bitter_rsids_with_ensembl

# Импорты: стандартные, сторонние, локальные
import os
import streamlit as st
from app.parsers import load_genotype_from_txt, load_genotype_from_vcf

def run_app():
    st.title("Распознавание восприятия горького вкуса по генотипу")


    st.write("""
    **Инструкция:**
    1. Перетащите или выберите файл генотипа (txt или vcf).
    2. Если не видите нужный файл, убедитесь, что у него расширение `.txt` или `.vcf` и он находится в нужной папке.
    3. Примеры файлов можно найти в папке `app/data` вашего проекта.
    """)

    # Показываем список файлов-примеров
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    if os.path.exists(data_dir):
        st.markdown('**Примеры файлов для загрузки:**')
        for fname in os.listdir(data_dir):
            st.code(fname)

    uploaded_file = st.file_uploader("Загрузите файл генотипа (txt или vcf)", type=["txt", "vcf"])


    # Кэшируем результат аннотации SNP для повторного использования

    @st.cache_data(show_spinner=False)
    def get_bitter_snp_annotations_sync():
        loop = asyncio.new_event_loop()
        try:
            asyncio.set_event_loop(loop)
            result = loop.run_until_complete(_get_bitter_snp_annotations_async())
        finally:
            loop.close()
        return result

    async def _get_bitter_snp_annotations_async():
        async with httpx.AsyncClient(timeout=30) as client:
            df = await annotate_bitter_rsids_with_ensembl(client)
            return df

    def analyze_taste(genotype: dict) -> str:
        """
        Анализирует чувствительность к горькому вкусу по генотипу и аннотированным rsID.
        """
        df = get_bitter_snp_annotations_sync()
        if df.empty:
            return "Не удалось получить аннотированные SNP для анализа."
        found = []
        for _, row in df.iterrows():
            rsid = row['rsid']
            allele = genotype.get(rsid, '')
            found.append((rsid, allele))
        if not any(allele for _, allele in found):
            return (
                "Не удалось найти необходимые SNP (" +
                ", ".join(df['rsid']) + ") в файле."
            )
        # Упрощённая логика: если хотя бы один T в любом из SNP — чувствителен
        if any('T' in allele for _, allele in found):
            return "Вы чувствительны к горькому вкусу (TAS2R38)."
        else:
            return "Вы менее чувствительны к горькому вкусу (TAS2R38)."

    if st.button('Показать аннотированные ключевые SNP TAS2R38'):
        with st.spinner('Аннотируем через Ensembl...'):
            df = get_bitter_snp_annotations_sync()
            st.dataframe(df)
            if not df.empty:
                fig = px.bar(
                    df,
                    x="protein_position",
                    y="rsid",
                    title="Ключевые SNP TAS2R38 по позициям белка",
                    labels={"protein_position": "Позиция в белке", "rsid": "rsID"}
                )
                st.plotly_chart(fig, use_container_width=True)

    if uploaded_file is not None:
        filename = uploaded_file.name
        file_ext = os.path.splitext(filename)[1].lower()
        st.info(f"Загружен файл: {filename} (расширение: {file_ext})")
        # Сохраняем временно файл
        temp_path = f"temp_uploaded{file_ext}"
        with open(temp_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
        try:
            # Определяем тип файла и парсим
            if file_ext == ".txt":
                st.write("Парсим как TXT...")
                genotype = load_genotype_from_txt(temp_path)
            elif file_ext == ".vcf":
                st.write("Парсим как VCF...")
                genotype = load_genotype_from_vcf(temp_path)
            else:
                st.error("Неподдерживаемый формат файла.")
                genotype = None
            st.write(f"Результат парсинга: {genotype}")
            if genotype is not None:
                result = analyze_taste(genotype)
                st.success(result)
        except Exception as e:
            st.error(f"Ошибка при обработке файла: {e}")
        finally:
            # Удаляем временный файл
            os.remove(temp_path)