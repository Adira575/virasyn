from virasyn import database


def test_entrez_summary_builds_pubmed_result(monkeypatch):
    payload = {
        "result": {
            "uids": ["12345"],
            "12345": {
                "title": "A synthetic example about viral proteins",
                "pubdate": "2025 Jan",
                "fulljournalname": "Journal of Test Virology",
                "authors": [{"name": "Example A"}, {"name": "Example B"}],
                "articleids": [{"idtype": "pubmed", "value": "12345"},
                               {"idtype": "doi", "value": "10.1000/example"}],
            },
        }
    }
    monkeypatch.setattr(database, "_get_ncbi_json", lambda *args, **kwargs: payload)
    rows = database._entrez_summary("pubmed", ["12345"])
    assert len(rows) == 1
    assert rows[0]["id"] == "12345"
    assert rows[0]["doi"] == "10.1000/example"
    assert rows[0]["authors"] == "Example A, Example B"
    assert rows[0]["url"] == "https://pubmed.ncbi.nlm.nih.gov/12345/"


def test_entrez_summary_builds_nucleotide_result(monkeypatch):
    payload = {
        "result": {
            "uids": ["6789"],
            "6789": {
                "accessionversion": "NC_000001.1",
                "title": "Example viral genome",
                "organism": "Synthetic virus",
                "slen": 4200,
                "updatedate": "2025/01/01",
            },
        }
    }
    monkeypatch.setattr(database, "_get_ncbi_json", lambda *args, **kwargs: payload)
    rows = database._entrez_summary("nuccore", ["6789"])
    assert rows[0]["accession"] == "NC_000001.1"
    assert rows[0]["organism"] == "Synthetic virus"
    assert rows[0]["length"] == 4200
