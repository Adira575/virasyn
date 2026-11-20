from unittest.mock import Mock, patch

import pytest

from virasyn.alphafold import _safe_structure_url, download_alphafold_structure, lookup_alphafold_prediction


def test_rejects_invalid_accession():
    with pytest.raises(ValueError):
        lookup_alphafold_prediction('not an accession!')


@patch('virasyn.alphafold.requests.get')
def test_lookup_returns_prediction(mock_get):
    response = Mock()
    response.status_code = 200
    response.json.return_value = [{'uniprotAccession': 'P0DTC2', 'globalMetricValue': 92.1}]
    response.raise_for_status.return_value = None
    mock_get.return_value = response
    result = lookup_alphafold_prediction('p0dtc2')
    assert result[0]['uniprotAccession'] == 'P0DTC2'
    mock_get.assert_called_once()


@patch('virasyn.alphafold.requests.get')
def test_lookup_404_returns_empty_list(mock_get):
    response = Mock()
    response.status_code = 404
    mock_get.return_value = response
    assert lookup_alphafold_prediction('P0DTC2') == []


def test_rejects_untrusted_structure_url():
    with pytest.raises(ValueError):
        _safe_structure_url({'pdbUrl': 'https://example.com/model.pdb'}, 'PDB')


@patch('virasyn.alphafold.requests.get')
def test_downloads_pdb_coordinates(mock_get):
    response = Mock()
    response.text = 'HEADER    ALPHAFOLD MODEL\nATOM      1  CA  ALA A   1'
    response.raise_for_status.return_value = None
    mock_get.return_value = response
    text, filename = download_alphafold_structure(
        {'uniprotAccession': 'P0DTC2', 'pdbUrl': 'https://alphafold.ebi.ac.uk/files/AF-P0DTC2-F1-model_v4.pdb'},
        'PDB',
    )
    assert 'ATOM' in text
    assert filename == 'alphafold_P0DTC2.pdb'
