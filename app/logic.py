
import asyncio
from typing import List, Dict, Any, Generator, Optional

import pandas as pd
import httpx

from app.constants import (
    BITTER_GENE, BITTER_POSITIONS, GENE_LOOKUP_URL, VARIATION_REGION_URL,
    VEP_ID_URL, ENSEMBL_HEADERS
)

async def get_gene_coordinates(gene_name: str, client: httpx.AsyncClient) -> Optional[Dict[str, Any]]:
    """
    Получить координаты гена по его названию через Ensembl REST API.

    :param gene_name: Название гена (например, 'TAS2R38').
    :type gene_name: str
    :return: Словарь с координатами гена или None при ошибке.
    :rtype: dict | None
    """
    gene_url = GENE_LOOKUP_URL.format(gene_name=gene_name)
    resp = await client.get(gene_url, headers=ENSEMBL_HEADERS)
    if resp.status_code != 200:
        print(f'Не удалось получить координаты для гена {gene_name}')
        return None
    return resp.json()

async def get_rsids_by_gene_ensembl(gene_name: str, client: httpx.AsyncClient) -> List[str]:
    """
    Получить список rsID для вариаций в регионе гена через Ensembl REST API.

    :param gene_name: Название гена.
    :type gene_name: str
    :return: Список rsID.
    :rtype: List[str]
    """
    gene_data = await get_gene_coordinates(gene_name, client)
    if not gene_data:
        return []
    region = (
        f"{gene_data['seq_region_name']}:{gene_data['start']}-{gene_data['end']}"
    )
    var_url = VARIATION_REGION_URL.format(region=region)
    var_resp = await client.get(var_url, headers=ENSEMBL_HEADERS)
    if var_resp.status_code != 200:
        print(f'Не удалось получить вариации для региона {region}')
        return []
    vars = var_resp.json()
    rsids = [v['id'] for v in vars if v['id'].startswith('rs')]
    return rsids

def chunked(lst: List[Any], n: int) -> Generator[List[Any], None, None]:
    """
    Разбить список на чанки по n элементов.

    :param lst: Исходный список.
    :type lst: List[Any]
    :param n: Размер чанка.
    :type n: int
    :yield: Чанк из n элементов.
    :rtype: Generator[List[Any], None, None]
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def filter_missense_tas2r38(entry: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Фильтрует аннотацию по missense_variant и нужным позициям для TAS2R38.

    :param entry: Аннотация варианта из VEP.
    :type entry: dict
    :return: Словарь с rsid, consequence_terms, protein_position или None, если не подходит.
    :rtype: dict | None
    """
    rsid = entry.get('id')
    for transcript in entry.get('transcript_consequences', []):
        if transcript.get('gene_symbol') == BITTER_GENE:
            cons = transcript.get('consequence_terms', [])
            pos = transcript.get('protein_start')
            if 'missense_variant' in cons and pos in BITTER_POSITIONS:
                return {
                    'rsid': rsid,
                    'consequence_terms': cons,
                    'protein_position': pos
                }
    return None

async def annotate_rsids_with_ensembl(rsids: List[str], client: httpx.AsyncClient, gene: str = BITTER_GENE) -> pd.DataFrame:
    """
    Аннотировать список rsID через Ensembl VEP API, фильтруя по missense_variant и позициям 49, 262, 296.

    :param rsids: Список rsID.
    :type rsids: List[str]
    :param gene: Название гена (по умолчанию TAS2R38).
    :type gene: str
    :return: DataFrame с rsid, consequence_terms, protein_position.
    :rtype: pd.DataFrame
    """
    results = []
    tasks = []
    for chunk in chunked(rsids, 200):
        tasks.append(
            client.post(
                VEP_ID_URL, headers=ENSEMBL_HEADERS, json={"ids": chunk}
            )
        )
    responses = await asyncio.gather(*tasks)
    for resp in responses:
        if resp.status_code != 200:
            print(f'Ошибка batch-запроса')
            continue
        data = resp.json()
        for entry in data:
            filtered = filter_missense_tas2r38(entry)
            if filtered:
                results.append(filtered)
    df = pd.DataFrame(results)
    df = df.drop_duplicates(subset=['protein_position'])
    return df

async def annotate_bitter_rsids_with_ensembl(client: httpx.AsyncClient) -> pd.DataFrame:
    """
    Получить и аннотировать ключевые SNP для TAS2R38 (горький вкус).

    :return: DataFrame с rsid, consequence_terms, protein_position.
    :rtype: pd.DataFrame
    """
    all_rsids = await get_rsids_by_gene_ensembl(BITTER_GENE, client)
    df = await annotate_rsids_with_ensembl(all_rsids, client)
    return df

async def main_async():
    async with httpx.AsyncClient(timeout=30) as client:
        df = await annotate_bitter_rsids_with_ensembl(client)
        print(df)

if __name__ == "__main__":
    asyncio.run(main_async())
