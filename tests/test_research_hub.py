from virasyn import database


def test_europe_pmc_record_mapping(monkeypatch):
    payload = {
        "hitCount": 1,
        "resultList": {"result": [{
            "title": "Phage receptor binding study",
            "firstPublicationDate": "2025-04-10",
            "journalTitle": "Virology",
            "authorString": "Example A, Example B",
            "pmid": "123",
            "pmcid": "PMC123",
            "doi": "10.1000/test",
            "isOpenAccess": "Y",
        }]},
    }
    monkeypatch.setattr(database, "_get_json", lambda *args, **kwargs: payload)
    result = database.search_europe_pmc("phage receptor", page=1, limit=10)
    assert result["count"] == 1
    record = result["results"][0]
    assert record["pmid"] == "123"
    assert record["doi"] == "10.1000/test"
    assert record["is_open_access"] is True


def test_openalex_record_mapping(monkeypatch):
    payload = {
        "meta": {"count": 1},
        "results": [{
            "id": "https://openalex.org/W1",
            "doi": "https://doi.org/10.1000/test",
            "title": "Viral host interaction",
            "publication_year": 2024,
            "authorships": [{"author": {"display_name": "Researcher One"}}],
            "primary_location": {"source": {"display_name": "Example Journal"}},
            "cited_by_count": 7,
            "open_access": {"is_oa": True},
        }],
    }
    monkeypatch.setattr(database, "_get_json", lambda *args, **kwargs: payload)
    result = database.search_openalex("viral host", page=1, limit=10)
    assert result["count"] == 1
    assert result["results"][0]["authors"] == "Researcher One"
    assert result["results"][0]["cited_by_count"] == 7


def test_interpro_entry_mapping(monkeypatch):
    payload = {
        "meta": {"count": 1},
        "results": [{
            "metadata": {
                "accession": "IPR000001",
                "name": {"name": "Example domain"},
                "description": [{"text": "A test protein domain"}],
                "type": "domain",
            }
        }],
    }
    monkeypatch.setattr(database, "_get_json", lambda *args, **kwargs: payload)
    result = database.search_interpro("example domain", page=1, limit=10)
    assert result["count"] == 1
    assert result["results"][0]["accession"] == "IPR000001"
    assert result["results"][0]["title"] == "Example domain"


