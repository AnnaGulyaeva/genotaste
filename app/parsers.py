def load_genotype_from_txt(filepath: str) -> dict:
    """Загрузка генотипа из файла 23andMe или MyHeritage (txt).
    Ожидается, что файл содержит строки с rsID и аллелями, разделённые табуляцией или пробелами.
    Возвращает словарь: {rsID: аллель}
    """
    genotype_data = {}
    with open(filepath, 'r', encoding='utf-8') as file_txt:
        for line in file_txt:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                rsid = parts[0]
                allele = parts[-1]
                genotype_data[rsid] = allele
    return genotype_data

def load_genotype_from_vcf(filepath: str) -> dict:
    """Загрузка генотипа из VCF файла.
    Ожидается, что файл содержит строки с rsID в третьей колонке и генотипом в 10-й (или последней) колонке.
    Возвращает словарь: {rsID: аллель}
    """
    genotype_data = {}
    with open(filepath, 'r', encoding='utf-8') as file_vcf:
        for line in file_vcf:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 10:
                rsid = parts[2]
                genotype_info = parts[9]
                # Обычно генотип в формате 0|1, 1|1 и т.д. (GT:...)
                gt = genotype_info.split(':')[0]
                # Преобразуем 0/1/2 в аллели, если нужно, иначе просто сохраняем как есть
                genotype_data[rsid] = gt
    return genotype_data