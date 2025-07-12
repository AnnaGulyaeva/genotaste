


# Импорты: стандартные, сторонние, локальные
from io import StringIO
from pathlib import Path
import pandas as pd
from app.logic import get_snp_annotations


def load_target_snps() -> dict:
    """Получение аннотаций SNP через функцию logic.get_snp_annotations."""
    results = get_snp_annotations(['rs713598', 'rs1726866', 'rs10246939'])
    print(results)
    print(results.columns.tolist())
    # Пример преобразования DataFrame в словарь (можно адаптировать под нужды)
    snp_dict = {}
    for _, row in results.iterrows():
        snp_dict[row['refsnp_id']] = {
            'gene': row.get('associated_gene', ''),
            'allele': row.get('allele', ''),
            'chr': row.get('chr_name', ''),
            'pos': row.get('chrom_start', '')
        }
    return snp_dict


target_snps = {
    "rs713598": "TAS2R38",
    "rs1726866": "TAS2R38",
    "rs10246939": "TAS2R38"
}

file_path = Path('C:/Users/guliaeva/PycharmProjects/open_AI/data/snp.vcf')


def load_vcf_to_dataframe() -> pd.DataFrame:
    """Загрузка VCF в pandas."""
    with open(file_path) as file_data:
        lines = [line for line in file_data if not line.startswith('##')]
    
    return pd.read_csv(
        StringIO('\n'.join(lines)),
        sep='\t'
    )


def interpret_taste_profile(genotipe_series: list) -> str:
    """Интерпретация вкусового профиля на основе генотипов."""
    sensitive_alleles = {"rs713598": "1", "rs1726866": "1", "rs10246939": "1"}
    score = 0

    for snp, gt in genotipe_series.items():
        alleles = gt.replace('|', '/').split('/')
        score += alleles.count(sensitive_alleles[snp])

    if score >= 5:
        return 'Вы очень чувствительны к горечи (возможно, не любите брокколи 🥦)'
    elif score >= 3:
        return 'Умеренная чувствительность к горечи'
    else:
        return ('Слабая чувствительность к горечи '
                '(брокколи вам, скорее всего, по вкусу!)')


if __name__ == '__main__':
    print(load_target_snps())

    df = load_vcf_to_dataframe()

    filtered = df[df['ID'].isin(target_snps.keys())].copy()

    filtered['GT'] = filtered['SAMPLE'].str.split(':').str[0]
    genotype_series = filtered.set_index('ID')['GT']

    print("Ваши генотипы вкусового гена TAS2R38:")
    for rsid, gt in genotype_series.items():
        print(f'    {rsid} ({target_snps[rsid]}): {gt}')
    
    print("\nИнтерпретация:")
    print(interpret_taste_profile(genotype_series))